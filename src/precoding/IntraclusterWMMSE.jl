##########################################################################
# IntraclusterWMMSE
#
# This is a version of the original WMMSE algorithm from
# Shi, Razaviyayn, Luo, He, "An iteratively weighted MMSE approach to
# distributed sum-utility maximization for a MIMO interfering broadcast
# channel", IEEE Trans. Signal Process., vol. 59, no. 9, 4331-4340, 2011,
# which has been modified to take into account the intracluster
# interference. This is done by maximizing a lower bound, where Jensen's
# inequality has been applied to the logdet(E) term, w.r.t. the expectation
# of out-of-cluster interference.
#
# Note that his precoding methods in parallel optimizes the cluster SDMA
# rate, as well as the network SDMA rate.

immutable IntraclusterWMMSEState
    # Receive filters (not necessarily MMSE filters)
    U_cluster_sdma::Array{Matrix{Complex128},1}
    U_network_sdma::Array{Matrix{Complex128},1}

    # Precoders
    V_cluster_sdma::Array{Matrix{Complex128},1}
    V_network_sdma::Array{Matrix{Complex128},1}

    # utility-optimal weights
    Z_cluster_sdma::Array{Diagonal{Float64},1}
    Z_network_sdma::Array{Diagonal{Float64},1}

    # logdet(inv(E_cluster_sdma_full[k])) is the cluster SDMA achievable rate
    # for MS k, i.e. with full CSI-R.
    E_cluster_sdma_full::Array{Diagonal{Float64},1}

    # logdet(inv(E_network_sdma_full[k])) is the network SDMA achievable rate
    # for MS k, i.e. with full CSI-R.
    # (That is, E_network_sdma_full is the MMSE matrix.)
    E_network_sdma_full::Array{Diagonal{Float64},1}

    # logdet(inv(E_network_sdma_full_partial[k])) is either a rate lower bound
    # (if the receiver is somehow aware of the existence of the intracluster
    # interference) or an achievable rate (if the receiver is oblivious of the
    # intracluster interference), cf. Lapidoth1996.
    E_network_sdma_partial::Array{Diagonal{Float64},1}
end

NaiveIntraclusterWMMSE(channel, network) =
    IntraclusterWMMSE(channel, network, network_sdma_robustness=false)
RobustIntraclusterWMMSE(channel, network) =
    IntraclusterWMMSE(channel, network, network_sdma_robustness=true)

function IntraclusterWMMSE(channel, network; network_sdma_robustness::Bool=true)
    K = get_num_MSs(network)
    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    ds = get_num_streams(network); max_d = maximum(ds)
    assignment = get_assignment(network)

    aux_params = get_aux_precoding_params(network)
    @defaultize_param! aux_params "IntraclusterWMMSE:bisection_Gamma_cond" 1e10
    @defaultize_param! aux_params "IntraclusterWMMSE:bisection_singular_Gamma_mu_lower_bound" 1e-14
    @defaultize_param! aux_params "IntraclusterWMMSE:bisection_max_iters" 5e1
    @defaultize_param! aux_params "IntraclusterWMMSE:bisection_tolerance" 1e-3

    prelogs_cluster_sdma = get_aux_network_param(network, "prelogs_cluster_sdma")
    prelogs_network_sdma = get_aux_network_param(network, "prelogs_network_sdma")

    state = IntraclusterWMMSEState(
        Array(Matrix{Complex128}, K),
        Array(Matrix{Complex128}, K),
        initial_precoders(channel, Ps, sigma2s, ds, assignment, aux_params),
        initial_precoders(channel, Ps, sigma2s, ds, assignment, aux_params),
        Array(Diagonal{Float64}, K),
        Array(Diagonal{Float64}, K),
        Array(Diagonal{Float64}, K),
        Array(Diagonal{Float64}, K),
        Array(Diagonal{Float64}, K)
    )
    objective = Float64[]
    utilities = Array(Float64, K, max_d, aux_params["max_iters"])
    weighted_logdet_rates_cluster_sdma_full = Array(Float64, K, max_d, aux_params["max_iters"])
    weighted_logdet_rates_network_sdma_full = Array(Float64, K, max_d, aux_params["max_iters"])
    weighted_logdet_rates_network_sdma_partial = Array(Float64, K, max_d, aux_params["max_iters"])

    iters = 0; conv_crit = Inf
    while iters < aux_params["max_iters"]
        update_MSs!(state, channel, Ps, sigma2s, ds, assignment, network_sdma_robustness)
        iters += 1

        # Results after this iteration
        utilities[:,:,iters] = calculate_utilities(state, prelogs_cluster_sdma, prelogs_network_sdma)
        push!(objective, sum(utilities[:,:,iters]))
        weighted_logdet_rates_cluster_sdma_full[:,:,iters], weighted_logdet_rates_network_sdma_full[:,:,iters], weighted_logdet_rates_network_sdma_partial[:,:,iters] =
            calculate_weighted_logdet_rates(state, prelogs_cluster_sdma, prelogs_network_sdma)

        # Check convergence
        if iters >= 2
            conv_crit = abs(objective[end] - objective[end-1])/abs(objective[end-1])
            if conv_crit < aux_params["stop_crit"] || objective[end] == 0
                Lumberjack.debug("IntraclusterWMMSE converged.",
                    @Compat.Dict(
                        :num_iters => iters,
                        :final_objective => objective[end],
                        :conv_crit => conv_crit,
                        :stop_crit => aux_params["stop_crit"],
                        :max_iters => aux_params["max_iters"])
                )
                break
            end
        end

        # Begin next iteration, unless the loop will end
        if iters < aux_params["max_iters"]
            update_BSs!(state, channel, Ps,
                prelogs_cluster_sdma, prelogs_network_sdma,
                assignment, aux_params, network_sdma_robustness)
        end
    end
    if iters == aux_params["max_iters"]
        Lumberjack.debug("IntraclusterWMMSE did NOT converge.",
            @Compat.Dict(
                :num_iters => iters,
                :final_objective => objective[end],
                :conv_crit => conv_crit,
                :stop_crit => aux_params["stop_crit"],
                :max_iters => aux_params["max_iters"])
        )
    end

    results = PrecodingResults()
    if aux_params["output_protocol"] == :all_iterations
        results["objective"] = objective
        results["utilities"] = utilities
        results["weighted_logdet_rates_cluster_sdma_full"]    = weighted_logdet_rates_cluster_sdma_full
        results["weighted_logdet_rates_network_sdma_full"]    = weighted_logdet_rates_network_sdma_full
        results["weighted_logdet_rates_network_sdma_partial"] = weighted_logdet_rates_network_sdma_partial
    elseif aux_params["output_protocol"] == :final_iteration
        results["objective"] = objective[iters]
        results["utilities"] = utilities[:,:,iters]
        results["weighted_logdet_rates_cluster_sdma_full"]    = weighted_logdet_rates_cluster_sdma_full[:,:,iters]
        results["weighted_logdet_rates_network_sdma_full"]    = weighted_logdet_rates_network_sdma_full[:,:,iters]
        results["weighted_logdet_rates_network_sdma_partial"] = weighted_logdet_rates_network_sdma_partial[:,:,iters]
    end
    results["weighted_logdet_rates_full"]    = results["weighted_logdet_rates_cluster_sdma_full"] + results["weighted_logdet_rates_network_sdma_full"]
    results["weighted_logdet_rates_partial"] = results["weighted_logdet_rates_cluster_sdma_full"] + results["weighted_logdet_rates_network_sdma_partial"]
    return results
end

# Note that cluster SDMA and naive network SDMA are not equivalent, because
# they might have different prelog factors between users.
function update_MSs!(state::IntraclusterWMMSEState, channel::SinglecarrierChannel,
    Ps, sigma2s, ds, assignment, network_sdma_robustness)

    for i = 1:channel.I; for k in served_MS_ids(i, assignment)
        coordinators = coordinated_BS_ids(k, assignment)

        # Covariances
        Phi_cluster_sdma_full           = complex(sigma2s[k]*eye(channel.Ns[k])) # intercluster instantaneous CSI-R
        Phi_network_sdma_full           = complex(sigma2s[k]*eye(channel.Ns[k])) # intercluster and intracluster instantaneous CSI-R
        Phi_network_sdma_partial_naive  = complex(sigma2s[k]*eye(channel.Ns[k])) # intercluster instantaneous CSI-R; intracluster CSI-R ignored
        Phi_network_sdma_partial_robust = complex(sigma2s[k]*eye(channel.Ns[k])) # intercluster instantaneous CSI-R; intracluster statistical CSI-R

        # Intracluster CSI-R
        for j in coordinators; for l in served_MS_ids(j, assignment)
            F_cluster_sdma = channel.H[k,j]*state.V_cluster_sdma[l]
            FFh_cluster_sdma = F_cluster_sdma*F_cluster_sdma'
            Phi_cluster_sdma_full           += FFh_cluster_sdma

            F_network_sdma = channel.H[k,j]*state.V_network_sdma[l]
            FFh_network_sdma = F_network_sdma*F_network_sdma'
            Phi_network_sdma_full           += FFh_network_sdma
            Phi_network_sdma_partial_naive  += FFh_network_sdma
            Phi_network_sdma_partial_robust += FFh_network_sdma
        end; end

        # Intercluster CSI-R
        for j in setdiff(IntSet(1:channel.I), coordinators); for l in served_MS_ids(j, assignment)
            F_network_sdma = channel.H[k,j]*state.V_network_sdma[l]
            FFh_network_sdma = F_network_sdma*F_network_sdma'
            Phi_network_sdma_full += FFh_network_sdma

            FFhrob_network_sdma = channel.large_scale_fading_factor[k,j]*channel.large_scale_fading_factor[k,j]*trace(state.V_network_sdma[l]'*state.V_network_sdma[l])*eye(channel.Ns[k])
            Phi_network_sdma_partial_robust += FFhrob_network_sdma
        end; end

        ### CLUSTER SDMA ###

        # Desired channel
        F_cluster_sdma = channel.H[k,i]*state.V_cluster_sdma[k]

        # Optimal (MMSE) receiver
        state.U_cluster_sdma[k] = inv(eigfact(hermitianize!(Phi_cluster_sdma_full)))*F_cluster_sdma

        # MMSE and optimal weight
        state.E_cluster_sdma_full[k] = Diagonal(min(1, abs(diag((eye(ds[k]) - state.U_cluster_sdma[k]'*F_cluster_sdma)))))
        state.Z_cluster_sdma[k] = inv(state.E_cluster_sdma_full[k])

        ### NETWORK SDMA ###

        # Desired channel
        F_network_sdma = channel.H[k,i]*state.V_network_sdma[k]

        # Receivers
        U_network_sdma_full   = inv(eigfact(hermitianize!(Phi_network_sdma_full)))*F_network_sdma
        U_network_sdma_naive  = inv(eigfact(hermitianize!(Phi_network_sdma_partial_naive)))*F_network_sdma
        U_network_sdma_robust = inv(eigfact(hermitianize!(Phi_network_sdma_partial_robust)))*F_network_sdma

        # MSEs and optimal weight
        state.E_network_sdma_full[k]    = Diagonal(min(1, abs(diag(eye(ds[k]) - U_network_sdma_full'*F_network_sdma))))
        state.E_network_sdma_partial[k] = Diagonal(min(1, abs(diag(eye(ds[k]) - U_network_sdma_robust'*F_network_sdma))))
        if network_sdma_robustness
            state.U_network_sdma[k] = U_network_sdma_robust
            state.Z_network_sdma[k] = inv(state.E_network_sdma_partial[k])
        else
            state.U_network_sdma[k] = U_network_sdma_naive
            state.Z_network_sdma[k] = inv(Diagonal(min(1, abs(diag((eye(ds[k]) - U_network_sdma_naive'*F_network_sdma))))))
        end
    end; end
end

function update_BSs!(state::IntraclusterWMMSEState, channel::SinglecarrierChannel,
    Ps, prelogs_cluster_sdma, prelogs_network_sdma,
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
                GZGh_cluster_sdma = prelogs_cluster_sdma[l]*G_cluster_sdma*state.Z_cluster_sdma[l]*G_cluster_sdma'
                Gamma_cluster_sdma_full           += GZGh_cluster_sdma

                G_network_sdma = channel.H[l,i]'*state.U_network_sdma[l]
                GZGh_network_sdma = prelogs_network_sdma[l]*G_network_sdma*state.Z_network_sdma[l]*G_network_sdma'
                Gamma_network_sdma_partial_naive  += GZGh_network_sdma
                Gamma_network_sdma_partial_robust += GZGh_network_sdma
            else
                # Intracluster CSI-T
                GZGhrob_network_sdma = prelogs_network_sdma[l]*channel.large_scale_fading_factor[l,i]*channel.large_scale_fading_factor[l,i]*trace(state.U_network_sdma[l]'*state.U_network_sdma[l]*state.Z_network_sdma[l])*eye(channel.Ms[i])
                Gamma_network_sdma_partial_robust += GZGhrob_network_sdma
            end
        end; end

        ### CLUSTER SDMA ###

        # Find optimal Lagrange multiplier
        mu_cluster_sdma, Gamma_eigen_cluster_sdma =
            optimal_mu(i, Gamma_cluster_sdma_full,               state.U_cluster_sdma, state.Z_cluster_sdma, prelogs_cluster_sdma, channel, Ps, assignment, aux_params)

        # Precoders (reuse EVD)
        for k in served_MS_ids(i, assignment)
            state.V_cluster_sdma[k] =
                Gamma_eigen_cluster_sdma.vectors*Diagonal(1./(abs(Gamma_eigen_cluster_sdma.values) .+ mu_cluster_sdma))*Gamma_eigen_cluster_sdma.vectors'*channel.H[k,i]'*state.U_cluster_sdma[k]*state.Z_cluster_sdma[k]*prelogs_cluster_sdma[k]
        end

        ### NETWORK SDMA ###

        # Find optimal Lagrange multiplier
        if network_sdma_robustness
            mu_network_sdma, Gamma_eigen_network_sdma =
                optimal_mu(i, Gamma_network_sdma_partial_robust, state.U_network_sdma, state.Z_network_sdma, prelogs_network_sdma, channel, Ps, assignment, aux_params)
        else
            mu_network_sdma, Gamma_eigen_network_sdma =
                optimal_mu(i, Gamma_network_sdma_partial_naive,  state.U_network_sdma, state.Z_network_sdma, prelogs_network_sdma, channel, Ps, assignment, aux_params)
        end

        # Precoders (reuse EVD)
        for k in served_MS_ids(i, assignment)
            state.V_network_sdma[k] =
                Gamma_eigen_network_sdma.vectors*Diagonal(1./(abs(Gamma_eigen_network_sdma.values) .+ mu_network_sdma))*Gamma_eigen_network_sdma.vectors'*channel.H[k,i]'*state.U_network_sdma[k]*state.Z_network_sdma[k]*prelogs_network_sdma[k]
        end
    end
end

# Calculates mu based on bisection.
function optimal_mu(i, Gamma, U, Z, prelogs, channel::SinglecarrierChannel, Ps, assignment, aux_params)
    # Build stuff needed for bisector function
    bis_M = complex(zeros(channel.Ms[i], channel.Ms[i]))
    for k in served_MS_ids(i, assignment)
        F = prelogs[k]*channel.H[k,i]'*U[k]*Z[k]
        FFh = F*F'
        bis_M += FFh
    end
    Gamma_eigen = eigfact(hermitianize!(Gamma)); Gamma_eigen_values = abs(Gamma_eigen.values)
    bis_JMJ_diag = abs(diag(Gamma_eigen.vectors'*bis_M*Gamma_eigen.vectors))
    bis_length = channel.Ms[i]

    # mu lower bound
    if maximum(Gamma_eigen_values)/minimum(Gamma_eigen_values) < aux_params["IntraclusterWMMSE:bisection_Gamma_cond"]
        # Gamma is invertible
        mu_lower = 0.
    else
        mu_lower = aux_params["IntraclusterWMMSE:bisection_singular_Gamma_mu_lower_bound"]
    end

    if bisector_IntraclusterWMMSE(mu_lower, bis_length, bis_JMJ_diag, Gamma_eigen_values) <= Ps[i]
        # No bisection needed
        return mu_lower, Gamma_eigen
    else
        # mu upper bound
        mu_upper = sqrt(channel.Ms[i]/Ps[i]*maximum(bis_JMJ_diag)) - minimum(Gamma_eigen_values)
        if bisector_IntraclusterWMMSE(mu_upper, bis_length, bis_JMJ_diag, Gamma_eigen_values) > Ps[i]
            Lumberjack.error("Power bisection: infeasible mu upper bound.")
        end

        num_iters = 0
        while num_iters < aux_params["IntraclusterWMMSE:bisection_max_iters"]
            conv_crit = (Ps[i] - bisector_IntraclusterWMMSE(mu_upper, bis_length, bis_JMJ_diag, Gamma_eigen_values))/Ps[i]

            if conv_crit < aux_params["IntraclusterWMMSE:bisection_tolerance"]
                break
            else
                mu = (1/2)*(mu_lower + mu_upper)

                if bisector_IntraclusterWMMSE(mu, bis_length, bis_JMJ_diag, Gamma_eigen_values) < Ps[i]
                    # New point feasible, replace upper point
                    mu_upper = mu
                else
                    # New point not feasible, replace lower point
                    mu_lower = mu
                end
            end

            num_iters += 1
        end

        if num_iters == aux_params["IntraclusterWMMSE:bisection_max_iters"]
            Lumberjack.warn("Power bisection: reached max iterations.")
        end

        # The upper point is always feasible, therefore we use it
        return mu_upper, Gamma_eigen
    end
end

# Non-anonymous bisector function, for speedup.
function bisector_IntraclusterWMMSE(mu, bis_length, bis_JMJ_diag, Gamma_eigen_values)
    output = 0.
    for n = 1:bis_length
        den_sqrt = Gamma_eigen_values[n] + mu
        output += bis_JMJ_diag[n]/(den_sqrt*den_sqrt)
    end
    return output
end

function calculate_utilities(state::IntraclusterWMMSEState, prelogs_cluster_sdma, prelogs_network_sdma)
    K = length(state.Z_cluster_sdma)
    ds = Int[ size(state.Z_cluster_sdma[k], 2) for k = 1:K ]; max_d = maximum(ds)

    utilities = zeros(Float64, K, max_d)
    for k = 1:K; if ds[k] > 0
        utilities[k,1:ds[k]] =
            prelogs_cluster_sdma[k]*log2(abs(eigvals(state.Z_cluster_sdma[k]))) +
            prelogs_network_sdma[k]*log2(abs(eigvals(state.Z_network_sdma[k])))
    end; end

    return utilities
end
