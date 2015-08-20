##########################################################################
# IntraclusterLeakageMinimization
#
# This is a version of the original leakage minimization algorithm from
# Gomadam, Cadambe, Jafar, "A Distributed Numerical Approach to
# Interference Alignment and Applications to Wireless Interference Networks,"
# IEEE Trans. IT, vol. 57, no. 6, pp. 3309–3322, 2011,
# which has been modified to be unaware of the spatial characteristics
# of the out-of-cluster interference. It has also been modified to account
# for the overhead pre-log factors when the interference-free subspaces
# are calculated. This will still lead to IA solution "convergence", but
# does affect the rate achieved.

immutable IntraclusterLeakageMinimizationState
    U_cluster_sdma::Array{Matrix{Complex128},1}
    U_network_sdma::Array{Matrix{Complex128},1}

    V_cluster_sdma::Array{Matrix{Complex128},1}
    V_network_sdma::Array{Matrix{Complex128},1}

    E_cluster_sdma_full::Array{Diagonal{Float64},1}
    E_network_sdma_full::Array{Diagonal{Float64},1}
    E_network_sdma_partial::Array{Diagonal{Float64},1}
end

NaiveIntraclusterLeakageMinimization(channel, network) =
    IntraclusterLeakageMinimization(channel, network, network_sdma_robustness=false)
RobustIntraclusterLeakageMinimization(channel, network) =
    IntraclusterLeakageMinimization(channel, network, network_sdma_robustness=true)

function IntraclusterLeakageMinimization(channel, network; network_sdma_robustness::Bool=true)
    assignment = get_assignment(network)

    K = get_no_MSs(network)
    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    ds = get_no_streams(network); max_d = maximum(ds)
    alphas = get_user_priorities(network)
    aux_params = get_aux_precoding_params(network)

    prelogs_cluster_sdma = get_aux_network_param(network, "prelogs_cluster_sdma")
    prelogs_network_sdma = get_aux_network_param(network, "prelogs_network_sdma")

    state = IntraclusterLeakageMinimizationState(
        Array(Matrix{Complex128}, K),
        Array(Matrix{Complex128}, K),
        initial_precoders(channel, Ps, sigma2s, ds, assignment, aux_params),
        initial_precoders(channel, Ps, sigma2s, ds, assignment, aux_params),
        Array(Diagonal{Float64}, K),
        Array(Diagonal{Float64}, K),
        Array(Diagonal{Float64}, K)
    )
    objective = Float64[]
    weighted_logdet_rates_cluster_sdma_full = Array(Float64, K, max_d, aux_params["max_iters"])
    weighted_logdet_rates_network_sdma_full = Array(Float64, K, max_d, aux_params["max_iters"])
    weighted_logdet_rates_network_sdma_partial = Array(Float64, K, max_d, aux_params["max_iters"])

    iters = 0; conv_crit = Inf
    while iters < aux_params["max_iters"]
        update_MSs!(state, channel, Ps, sigma2s, ds, assignment, network_sdma_robustness)
        iters += 1

        # Results after this iteration
        weighted_logdet_rates_cluster_sdma_full[:,:,iters], weighted_logdet_rates_network_sdma_full[:,:,iters], weighted_logdet_rates_network_sdma_partial[:,:,iters] =
            calculate_weighted_logdet_rates(state, prelogs_cluster_sdma, prelogs_network_sdma)
        push!(objective,
            sum(Diagonal(prelogs_cluster_sdma)*weighted_logdet_rates_cluster_sdma_full[:,:,iters]) +
            sum(Diagonal(prelogs_network_sdma)*weighted_logdet_rates_network_sdma_full[:,:,iters]))

        # Check convergence
        if iters >= 2
            conv_crit = abs(objective[end] - objective[end-1])/abs(objective[end-1])
            if conv_crit < aux_params["stop_crit"] || objective[end] == 0
                Lumberjack.debug("IntraclusterLeakageMinimization converged.",
                    [ :num_iters => iters,
                      :final_objective => objective[end],
                      :conv_crit => conv_crit,
                      :stop_crit => aux_params["stop_crit"],
                      :max_iters => aux_params["max_iters"] ]
                )
                break
            end
        end

        # Begin next iteration, unless the loop will end
        if iters < aux_params["max_iters"]
            update_BSs!(state, channel, Ps, ds,
                prelogs_cluster_sdma, prelogs_network_sdma,
                assignment, aux_params, network_sdma_robustness)
        end
    end
    if iters == aux_params["max_iters"]
        Lumberjack.debug("IntraclusterLeakageMinimization did NOT converge.",
            [ :num_iters => iters,
              :final_objective => objective[end],
              :conv_crit => conv_crit,
              :stop_crit => aux_params["stop_crit"],
              :max_iters => aux_params["max_iters"] ]
        )
    end

    results = PrecodingResults()
    if aux_params["output_protocol"] == :all_iterations
        results["objective"] = objective
        results["weighted_logdet_rates_cluster_sdma_full"]    = weighted_logdet_rates_cluster_sdma_full
        results["weighted_logdet_rates_network_sdma_full"]    = weighted_logdet_rates_network_sdma_full
        results["weighted_logdet_rates_network_sdma_partial"] = weighted_logdet_rates_network_sdma_partial
    elseif aux_params["output_protocol"] == :final_iteration
        results["objective"] = objective[iters]
        results["weighted_logdet_rates_cluster_sdma_full"]    = weighted_logdet_rates_cluster_sdma_full[:,:,iters]
        results["weighted_logdet_rates_network_sdma_full"]    = weighted_logdet_rates_network_sdma_full[:,:,iters]
        results["weighted_logdet_rates_network_sdma_partial"] = weighted_logdet_rates_network_sdma_partial[:,:,iters]
    end
    results["weighted_logdet_rates_full"]    = results["weighted_logdet_rates_cluster_sdma_full"] + results["weighted_logdet_rates_network_sdma_full"]
    results["weighted_logdet_rates_partial"] = results["weighted_logdet_rates_cluster_sdma_full"] + results["weighted_logdet_rates_network_sdma_partial"]
    return results
end

function update_MSs!(state::IntraclusterLeakageMinimizationState, channel::SinglecarrierChannel,
    Ps, sigma2s, ds, assignment, network_sdma_robustness)

    for i = 1:channel.I; for k in served_MS_ids(i, assignment)
        coordinators = coordinated_BS_ids(k, assignment)

        # Covariances
        Psi_cluster_sdma_full           = complex(zeros(channel.Ns[k], channel.Ns[k])) # intercluster instantaneous CSI-R
        Psi_network_sdma_full           = complex(zeros(channel.Ns[k], channel.Ns[k])) # intercluster and intracluster instantaneous CSI-R
        Psi_network_sdma_partial_naive  = complex(zeros(channel.Ns[k], channel.Ns[k])) # intercluster instantaneous CSI-R; intracluster CSI-R ignored
        Psi_network_sdma_partial_robust = complex(zeros(channel.Ns[k], channel.Ns[k])) # intercluster instantaneous CSI-R; intracluster statistical CSI-R

        # Intracluster CSI-R
        for j in coordinators; for l in served_MS_ids(j, assignment)
            if (j, l) == (i, k); continue; end

            F_cluster_sdma = channel.H[k,j]*state.V_cluster_sdma[l]
            FFh_cluster_sdma = F_cluster_sdma*F_cluster_sdma'
            Psi_cluster_sdma_full           += FFh_cluster_sdma

            F_network_sdma = channel.H[k,j]*state.V_network_sdma[l]
            FFh_network_sdma = F_network_sdma*F_network_sdma'
            Psi_network_sdma_full           += FFh_network_sdma
            Psi_network_sdma_partial_naive  += FFh_network_sdma
            Psi_network_sdma_partial_robust += FFh_network_sdma
        end; end

        # Intercluster CSI-R
        for j in setdiff(IntSet(1:channel.I), coordinators); for l in served_MS_ids(j, assignment)
            F_network_sdma = channel.H[k,j]*state.V_network_sdma[l]
            FFh_network_sdma = F_network_sdma*F_network_sdma'
            Psi_network_sdma_full += FFh_network_sdma

            FFhrob_network_sdma = Ps[j]*channel.large_scale_fading_factor[k,j]*channel.large_scale_fading_factor[k,j]*eye(channel.Ns[k])
            Psi_network_sdma_partial_robust += FFhrob_network_sdma
        end; end

        ### CLUSTER SDMA ###

        # Desired channel
        F_cluster_sdma = channel.H[k,i]*state.V_cluster_sdma[k]

        # Optimal (MMSE) receiver and MSE
        Ummse_cluster_sdma = inv(eigfact(Hermitian(Psi_cluster_sdma_full + F_cluster_sdma*F_cluster_sdma')))*F_cluster_sdma
        state.E_cluster_sdma_full[k] = Diagonal(min(1, abs(diag((eye(ds[k]) - Ummse_cluster_sdma'*F_cluster_sdma)))))

        # Actual (zero-forcing) receiver
        state.U_cluster_sdma[k] = eigfact(Hermitian(Psi_cluster_sdma_full), 1:ds[k]).vectors

        ### NETWORK SDMA ###

        # Desired channel
        F_network_sdma = channel.H[k,i]*state.V_network_sdma[k]

        # Optimal (MMSE) receivers and MSEs
        Ummse_network_sdma_full = inv(eigfact(Hermitian(Psi_network_sdma_full + F_network_sdma*F_network_sdma')))*F_network_sdma
        state.E_network_sdma_full[k] = Diagonal(min(1, abs(diag(eye(ds[k]) - Ummse_network_sdma_full'*F_network_sdma))))
        Ummse_network_sdma_robust = inv(eigfact(Hermitian(Psi_network_sdma_partial_robust + F_network_sdma*F_network_sdma')))*F_network_sdma
        state.E_network_sdma_partial[k] = Diagonal(min(1, abs(diag(eye(ds[k]) - Ummse_network_sdma_robust'*F_network_sdma))))

        # Actual (zero-forcing) receivers
        if network_sdma_robustness
            state.U_network_sdma[k] = eigfact(Hermitian(Psi_network_sdma_partial_robust), 1:ds[k]).vectors
        else
            state.U_network_sdma[k] = eigfact(Hermitian(Psi_network_sdma_partial_naive), 1:ds[k]).vectors
        end
    end; end
end

function update_BSs!(state::IntraclusterLeakageMinimizationState,
    channel::SinglecarrierChannel, Ps, ds,
    prelogs_cluster_sdma, prelogs_network_sdma,
    assignment, aux_params, network_sdma_robustness)

    for i in active_BSs(assignment)
        coordinatees = coordinated_MS_ids(i, assignment)

        # Covariances
        Gamma_cluster_sdma_full           = complex(zeros(channel.Ms[i], channel.Ms[i]))
        Gamma_network_sdma_partial_naive  = complex(zeros(channel.Ms[i], channel.Ms[i]))
        Gamma_network_sdma_partial_robust = complex(zeros(channel.Ms[i], channel.Ms[i]))

        for j = 1:channel.I; for l in served_MS_ids(j, assignment)
            if l in coordinatees
                # Intercluster CSI-T
                G_cluster_sdma = channel.H[l,i]'*state.U_cluster_sdma[l]
                GGh_cluster_sdma = prelogs_cluster_sdma[l]*G_cluster_sdma*G_cluster_sdma'
                Gamma_cluster_sdma_full           += GGh_cluster_sdma

                G_network_sdma = channel.H[l,i]'*state.U_network_sdma[l]
                GGh_network_sdma = prelogs_network_sdma[l]*G_network_sdma*G_network_sdma'
                Gamma_network_sdma_partial_naive  += GGh_network_sdma
                Gamma_network_sdma_partial_robust += GGh_network_sdma
            else
                # Intracluster CSI-T
                GGhrob_network_sdma = Ps[i]*prelogs_network_sdma[l]*channel.large_scale_fading_factor[l,i]*channel.large_scale_fading_factor[l,i]*eye(channel.Ms[i])
                Gamma_network_sdma_partial_robust += GGhrob_network_sdma
            end
        end; end

        ### CLUSTER SDMA ###

        served = served_MS_ids(i, assignment)
        Nserved = length(served)
        for k in served
            # Desired channel
            G_cluster_sdma = channel.H[k,i]'*state.U_cluster_sdma[k]

            state.V_cluster_sdma[k] = sqrt(Ps[i]/(Nserved*ds[k]))*eigfact(Hermitian(Gamma_cluster_sdma_full - G_cluster_sdma*G_cluster_sdma'), 1:ds[k]).vectors
        end

        ### NETWORK SDMA ###

        served = served_MS_ids(i, assignment)
        Nserved = length(served)
        for k in served
            # Desired channel
            G_network_sdma = channel.H[k,i]'*state.U_network_sdma[k]

            if network_sdma_robustness
                state.V_network_sdma[k] = sqrt(Ps[i]/(Nserved*ds[k]))*eigfact(Hermitian(Gamma_network_sdma_partial_robust - G_network_sdma*G_network_sdma'), 1:ds[k]).vectors
            else
                state.V_network_sdma[k] = sqrt(Ps[i]/(Nserved*ds[k]))*eigfact(Hermitian(Gamma_network_sdma_partial_naive - G_network_sdma*G_network_sdma'), 1:ds[k]).vectors
            end
        end
    end
end
