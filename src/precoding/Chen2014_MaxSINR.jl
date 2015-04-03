##########################################################################
# Chen2014_MaxSINR
#
# Out-of-cluster interference robustified version of MaxSINR from
# Chen, Cheng, "Clustering for Interference Alignment in Multiuser
# Interference Network", IEEE Trans. VT, vol. 63, no. 6, pp. 2613-2624, 2014.

immutable Chen2014_MaxSINRState
    U::Array{Matrix{Complex128},1}
    W::Array{Hermitian{Complex128},1} # these are only used for rate calculations
    V::Array{Matrix{Complex128},1}
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

    state = Chen2014_MaxSINRState(
        initial_receivers(channel, Ps, sigma2s, ds, assignment, aux_params),
        Array(Hermitian{Complex128}, K),
        initial_precoders(channel, Ps, sigma2s, ds, assignment, aux_params)
    )
    objective = Float64[]
    logdet_rates = Array(Float64, K, max_d, aux_params["max_iters"])
    MMSE_rates = Array(Float64, K, max_d, aux_params["max_iters"])
    weighted_logdet_rates = Array(Float64, K, max_d, aux_params["max_iters"])
    weighted_MMSE_rates = Array(Float64, K, max_d, aux_params["max_iters"])
    allocated_power = Array(Float64, K, max_d, aux_params["max_iters"])

    iters = 0; conv_crit = Inf
    while iters < aux_params["max_iters"]
        update_MSs!(state, channel, Ps, sigma2s, assignment, robustness)
        iters += 1

        # Results after this iteration
        logdet_rates[:,:,iters] = calculate_logdet_rates(state)
        MMSE_rates[:,:,iters] = calculate_MMSE_rates(state)
        weighted_logdet_rates[:,:,iters] = calculate_weighted_logdet_rates(state, alphas)
        push!(objective, sum(weighted_logdet_rates[:,:,iters]))
        weighted_MMSE_rates[:,:,iters] = calculate_weighted_MMSE_rates(state, alphas)
        allocated_power[:,:,iters] = calculate_allocated_power(state)

        # Check convergence
        if iters >= 2
            conv_crit = abs(objective[end] - objective[end-1])/abs(objective[end-1])
            if conv_crit < aux_params["stop_crit"]
                Lumberjack.debug("Chen2014_MaxSINR converged.",
                    [ :no_iters => iters,
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
            update_BSs!(state, channel, Ps, sigma2s, assignment, aux_params, robustness)
        end
    end
    if iters == aux_params["max_iters"]
        Lumberjack.debug("Chen2014_MaxSINR did NOT converge.",
            [ :no_iters => iters,
              :final_objective => objective[end],
              :conv_crit => conv_crit,
              :stop_crit => aux_params["stop_crit"],
              :max_iters => aux_params["max_iters"] ]
        )
    end

    results = PrecodingResults()
    if aux_params["output_protocol"] == :all_iterations
        results["objective"] = objective
        results["logdet_rates"] = logdet_rates
        results["MMSE_rates"] = MMSE_rates
        results["weighted_logdet_rates"] = weighted_logdet_rates
        results["weighted_MMSE_rates"] = weighted_MMSE_rates
        results["allocated_power"] = allocated_power
    elseif aux_params["output_protocol"] == :final_iteration
        results["objective"] = objective[iters]
        results["logdet_rates"] = logdet_rates[:,:,iters]
        results["MMSE_rates"] = MMSE_rates[:,:,iters]
        results["weighted_logdet_rates"] = weighted_logdet_rates[:,:,iters]
        results["weighted_MMSE_rates"] = weighted_MMSE_rates[:,:,iters]
        results["allocated_power"] = allocated_power[:,:,iters]
    end
    return results
end

function update_MSs!(state::Chen2014_MaxSINRState,
    channel::SinglecarrierChannel, Ps, sigma2s, assignment, robustness)

    ds = [ size(state.V[k], 2) for k = 1:channel.K ]

    for i = 1:channel.I; for k in served_MS_ids(i, assignment)
        coordinators = coordinated_BS_ids(k, assignment)

        # Covariances
        Phi_perfect = Hermitian(complex(sigma2s[k]*eye(channel.Ns[k])))
        Phi_imperfect = Hermitian(complex(sigma2s[k]*eye(channel.Ns[k])))
        for j in coordinators; for l in served_MS_ids(j, assignment)
            #Phi += Hermitian(channel.H[k,j]*(state.V[l]*state.V[l]')*channel.H[k,j]')
            Base.LinAlg.BLAS.herk!(Phi_perfect.uplo, 'N', complex(1.), channel.H[k,j]*state.V[l], complex(1.), Phi_perfect.S)
            Base.LinAlg.BLAS.herk!(Phi_imperfect.uplo, 'N', complex(1.), channel.H[k,j]*state.V[l], complex(1.), Phi_imperfect.S)
        end; end
        for j in setdiff(1:channel.I, coordinators); for l in served_MS_ids(j, assignment)
            Base.LinAlg.BLAS.herk!(Phi_perfect.uplo, 'N', complex(1.), channel.H[k,j]*state.V[l], complex(1.), Phi_perfect.S)
            robustness && (Phi_imperfect += Hermitian(complex(channel.large_scale_fading_factor[k,j]^2*Ps[j]*eye(channel.Ns[k]))))
        end; end

        # Per-stream receivers
        for n = 1:ds[k]
            Psi_plus_n = Hermitian(
                Base.LinAlg.BLAS.herk!(Phi_imperfect.uplo, 'N', complex(-1.), channel.H[k,i]*state.V[k][:,n], complex(1.), copy(Phi_imperfect.S)),
                Phi_imperfect.uplo)
            u = Psi_plus_n\channel.H[k,i]*state.V[k][:,n]
            state.U[k][:,n] = u/norm(u,2)
        end

        # True MSE weights (for rate calculation only)
        F = channel.H[k,i]*state.V[k]
        Ummse = Phi_perfect\F
        state.W[k] = Hermitian((eye(ds[k]) - Ummse'*F)\eye(ds[k]))
    end; end
end

function update_BSs!(state::Chen2014_MaxSINRState,
    channel::SinglecarrierChannel, Ps, sigma2s,
    assignment, aux_params, robustness)

    ds = [ size(state.V[k], 2) for k = 1:channel.K ]

    for i in active_BSs(assignment)
        coordinators = coordinated_MS_ids(i, assignment)

        # Covariance
        Gamma = Hermitian(complex(zeros(channel.Ms[i], channel.Ms[i])))
        for j = 1:channel.I; for l in served_MS_ids(j, assignment)
            if l in coordinators
                Gamma += Hermitian(channel.H[l,i]'*(state.U[l]*state.U[l]')*channel.H[l,i])
            else
                if robustness
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
