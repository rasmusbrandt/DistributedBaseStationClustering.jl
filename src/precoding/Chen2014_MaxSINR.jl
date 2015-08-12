##########################################################################
# Chen2014_MaxSINR
#
# Out-of-cluster interference robustified version of MaxSINR from
# Chen, Cheng, "Clustering for Interference Alignment in Multiuser
# Interference Network", IEEE Trans. VT, vol. 63, no. 6, pp. 2613-2624, 2014.

immutable Chen2014_MaxSINRState
    U::Array{Matrix{Complex128},1}
    V::Array{Matrix{Complex128},1}

    # See IntraclusterWMMSEState for description of these.
    E_full::Array{Diagonal{Float64},1}
    E_partial::Array{Diagonal{Float64},1}
end

NaiveChen2014_MaxSINR(channel, network) =
    Chen2014_MaxSINR(channel, network, robustness=false)
RobustChen2014_MaxSINR(channel, network) =
    Chen2014_MaxSINR(channel, network, robustness=true)

function Chen2014_MaxSINR(channel, network; robustness::Bool=true)
    assignment = get_assignment(network)

    K = get_no_MSs(network)
    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    ds = get_no_streams(network); max_d = maximum(ds)
    alphas = get_user_priorities(network)
    aux_params = get_aux_precoding_params(network)

    clustering_type = get_aux_assignment_param(network, "clustering_type")

    state = Chen2014_MaxSINRState(
        initial_receivers(channel, Ps, sigma2s, ds, assignment, aux_params),
        initial_precoders(channel, Ps, sigma2s, ds, assignment, aux_params),
        Array(Diagonal{Float64}, K),
        Array(Diagonal{Float64}, K)
    )
    objective = Float64[]
    weighted_logdet_rates_full = Array(Float64, K, max_d, aux_params["max_iters"])
    weighted_logdet_rates_partial = Array(Float64, K, max_d, aux_params["max_iters"])
    allocated_power = Array(Float64, K, max_d, aux_params["max_iters"])

    iters = 0; conv_crit = Inf
    while iters < aux_params["max_iters"]
        update_MSs!(state, channel, Ps, sigma2s, assignment,
            clustering_type, robustness)
        iters += 1

        # Results after this iteration
        weighted_logdet_rates_full[:,:,iters], weighted_logdet_rates_partial[:,:,iters] = calculate_weighted_logdet_rates(state, alphas)
        push!(objective, sum(weighted_logdet_rates_full[:,:,iters]))
        allocated_power[:,:,iters] = calculate_allocated_power(state)

        # Check convergence
        if iters >= 2
            conv_crit = abs(objective[end] - objective[end-1])/abs(objective[end-1])
            if conv_crit < aux_params["stop_crit"] || objective[end] == 0
                Lumberjack.debug("Chen2014_MaxSINR converged.",
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
            update_BSs!(state, channel, Ps, sigma2s, assignment, aux_params,
                clustering_type, robustness)
        end
    end
    if iters == aux_params["max_iters"]
        Lumberjack.debug("Chen2014_MaxSINR did NOT converge.",
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
        results["weighted_logdet_rates_full"] = weighted_logdet_rates_full
        results["weighted_logdet_rates_partial"] = weighted_logdet_rates_partial
        results["allocated_power"] = allocated_power
    elseif aux_params["output_protocol"] == :final_iteration
        results["objective"] = objective[iters]
        results["weighted_logdet_rates_full"] = weighted_logdet_rates_full[:,:,iters]
        results["weighted_logdet_rates_partial"] = weighted_logdet_rates_partial[:,:,iters]
        results["allocated_power"] = allocated_power[:,:,iters]
    end
    return results
end

function update_MSs!(state::Chen2014_MaxSINRState,
    channel::SinglecarrierChannel, Ps, sigma2s, assignment,
    clustering_type, robustness)

    ds = [ size(state.V[k], 2) for k = 1:channel.K ]

    for i = 1:channel.I; for k in served_MS_ids(i, assignment)
        coordinators = coordinated_BS_ids(k, assignment)

        # Covariances
        Phi_full           = Hermitian(complex(sigma2s[k]*eye(channel.Ns[k]))) # both intercluster and intracluster channels
        Phi_partial_naive  = Hermitian(complex(sigma2s[k]*eye(channel.Ns[k]))) # intercluster instantaneous CSI-R; intracluster ignored
        Phi_partial_robust = Hermitian(complex(sigma2s[k]*eye(channel.Ns[k]))) # intercluster instantaneous CSI-R; intracluster statistical CSI-R

        # Intercluster CSI-R
        for j in coordinators; for l in served_MS_ids(j, assignment)
            Fkjl = channel.H[k,j]*state.V[l]
            Base.LinAlg.BLAS.herk!(Phi_full.uplo, 'N', complex(1.), Fkjl, complex(1.), Phi_full.S)
            Base.LinAlg.BLAS.herk!(Phi_partial_naive.uplo, 'N', complex(1.), Fkjl, complex(1.), Phi_partial_naive.S)
            Base.LinAlg.BLAS.herk!(Phi_partial_robust.uplo, 'N', complex(1.), Fkjl, complex(1.), Phi_partial_robust.S)
        end; end

        # Intracluster CSI-R (only applicable if we are using spectrum_sharing!)
        if clustering_type == :spectrum_sharing
            for j in setdiff(IntSet(1:channel.I), coordinators); for l in served_MS_ids(j, assignment)
                Base.LinAlg.BLAS.herk!(Phi_full.uplo, 'N', complex(1.), channel.H[k,j]*state.V[l], complex(1.), Phi_full.S)
                Phi_partial_robust += Hermitian(complex(channel.large_scale_fading_factor[k,j]^2*Ps[j]*eye(channel.Ns[k])))
            end; end
        end

        # Desired channel
        Fiki = channel.H[k,i]*state.V[k]

        # Utility-optimal receiver (per-stream processing)
        if robustness
            for n = 1:ds[k]
                Psi_plus_n = Hermitian(
                    Base.LinAlg.BLAS.herk!(Phi_partial_robust.uplo, 'N', complex(-1.), channel.H[k,i]*state.V[k][:,n], complex(1.), copy(Phi_partial_robust.S)),
                    Phi_partial_robust.uplo)
                u = Psi_plus_n\channel.H[k,i]*state.V[k][:,n]
                state.U[k][:,n] = u/norm(u,2)
            end
        else
            for n = 1:ds[k]
                Psi_plus_n = Hermitian(
                    Base.LinAlg.BLAS.herk!(Phi_partial_naive.uplo, 'N', complex(-1.), channel.H[k,i]*state.V[k][:,n], complex(1.), copy(Phi_partial_naive.S)),
                    Phi_partial_naive.uplo)
                u = Psi_plus_n\channel.H[k,i]*state.V[k][:,n]
                state.U[k][:,n] = u/norm(u,2)
            end
        end

        # Full CSI, i.e. intracluster CSI is estimated in a final training stage.
        # This is an achievable rate.
        state.E_full[k] = Diagonal(min(1, abs(diag(eye(ds[k]) - (Phi_full\Fiki)'*Fiki))))

        # Partial CSI, i.e. intracluster CSI is _NOT_ estimated.
        # This is either a rate bound (if the receiver is aware that intracluster interference exists),
        # or an achievable rate (if the receiver is oblivious of the intracluster interference). Cf. Lapidoth1996.
        state.E_partial[k] = Diagonal(min(1, abs(diag(eye(ds[k]) - (Phi_partial_robust\Fiki)'*Fiki))))
    end; end
end

function update_BSs!(state::Chen2014_MaxSINRState,
    channel::SinglecarrierChannel, Ps, sigma2s, assignment, aux_params,
    clustering_type, robustness)

    ds = [ size(state.V[k], 2) for k = 1:channel.K ]

    for i in active_BSs(assignment)
        coordinators = coordinated_MS_ids(i, assignment)

        # Covariance
        Gamma = Hermitian(complex(zeros(channel.Ms[i], channel.Ms[i])))
        for j = 1:channel.I; for l in served_MS_ids(j, assignment)
            if l in coordinators
                Gamma += Hermitian(channel.H[l,i]'*(state.U[l]*state.U[l]')*channel.H[l,i])
            else
                # Take out-of-cluster interference into account if we are using
                # spectrum sharing, and we actually want to be robust against
                # this type of generated interference.
                if clustering_type == :spectrum_sharing && robustness
                    Gamma += Hermitian(complex(channel.large_scale_fading_factor[l,i]^2*eye(channel.Ms[i])))
                end
            end
        end; end

        # Per-stream precoders
        served = served_MS_ids(i, assignment)
        Nserved = length(served)
        for k in served
            for n = 1:ds[k]
                Gamma_i_plus_n = Hermitian(
                    Base.LinAlg.BLAS.herk!(Gamma.uplo, 'N', complex(-1.), channel.H[k,i]'*state.U[k][:,n], complex(1.), copy(Gamma.S)) + (sigma2s[k]/Ps[i])*eye(channel.Ms[i]),
                    Gamma.uplo)
                v = Gamma_i_plus_n\channel.H[k,i]'*state.U[k][:,n]
                state.V[k][:,n] = sqrt(Ps[i]/(Nserved*ds[k]))*v/norm(v,2)
            end
        end
    end
end
