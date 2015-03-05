immutable IntraclusterLeakageMinimizationState
    U::Array{Matrix{Complex128},1}
    W::Array{Hermitian{Complex128},1} # for rate calculations
    V::Array{Matrix{Complex128},1}
end

NaiveIntraclusterLeakageMinimization(channel, network) = IntraclusterLeakageMinimization(channel, network, robustness=false)
RobustIntraclusterLeakageMinimization(channel, network) = IntraclusterLeakageMinimization(channel, network, robustness=true)

function IntraclusterLeakageMinimization(channel, network; robustness::Bool=true)
    assignment = get_assignment(network)

    K = get_no_MSs(network)
    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    ds = get_no_streams(network)
    aux_params = get_aux_precoding_params(network)

    state = IntraclusterLeakageMinimizationState(
        Array(Matrix{Complex128}, K),
        unity_MSE_weights(ds),
        initial_precoders(channel, Ps, sigma2s, ds, assignment, aux_params))
    objective = Float64[]
    logdet_rates = Array(Float64, K, maximum(ds), aux_params["max_iters"])
    MMSE_rates = Array(Float64, K, maximum(ds), aux_params["max_iters"])
    allocated_power = Array(Float64, K, maximum(ds), aux_params["max_iters"])

    iters = 0; conv_crit = Inf
    while iters < aux_params["max_iters"]
        update_MSs!(state, channel, Ps, sigma2s, assignment, robustness)
        iters += 1

        # Results after this iteration
        logdet_rates[:,:,iters] = calculate_logdet_rates(state)
        push!(objective, sum(logdet_rates[:,:,iters]))
        MMSE_rates[:,:,iters] = calculate_MMSE_rates(state)
        allocated_power[:,:,iters] = calculate_allocated_power(state)

        # Check convergence
        if iters >= 2
            conv_crit = abs(objective[end] - objective[end-1])/abs(objective[end-1])
            if conv_crit < aux_params["stop_crit"]
                Lumberjack.debug("IntraclusterLeakageMinimization converged.",
                    [ :no_iters => iters, :final_objective => objective[end],
                      :conv_crit => conv_crit, :stop_crit => aux_params["stop_crit"],
                      :max_iters => aux_params["max_iters"] ])
                break
            end
        end

        # Begin next iteration, unless the loop will end
        if iters < aux_params["max_iters"]
            update_BSs!(state, channel, Ps, assignment, aux_params, robustness)
        end
    end
    if iters == aux_params["max_iters"]
        Lumberjack.debug("IntraclusterLeakageMinimization did NOT converge.",
            [ :no_iters => iters, :final_objective => objective[end],
              :conv_crit => conv_crit, :stop_crit => aux_params["stop_crit"],
              :max_iters => aux_params["max_iters"] ])
    end

    results = PrecodingResults()
    if aux_params["output_protocol"] == :all_iterations
        results["objective"] = objective
        results["logdet_rates"] = logdet_rates
        results["MMSE_rates"] = MMSE_rates
        results["allocated_power"] = allocated_power
    elseif aux_params["output_protocol"] == :final_iteration
        results["objective"] = objective[iters]
        results["logdet_rates"] = logdet_rates[:,:,iters]
        results["MMSE_rates"] = MMSE_rates[:,:,iters]
        results["allocated_power"] = allocated_power[:,:,iters]
    end
    return results
end

function update_MSs!(state::IntraclusterLeakageMinimizationState,
    channel::SinglecarrierChannel, Ps, sigma2s, assignment, robustness)

    ds = [ size(state.W[k], 1) for k = 1:channel.K ]

    for i = 1:channel.I; for k in served_MS_ids(i, assignment)
        coordinators = coordinated_BS_ids(k, assignment)

        # Covariances
        Phi_perfect = Hermitian(complex(sigma2s[k]*eye(channel.Ns[k])))
        Psi_imperfect = Hermitian(complex(sigma2s[k]*eye(channel.Ns[k])))
        for j in coordinators; for l in served_MS_ids(j, assignment)
            #Phi += Hermitian(channel.H[k,j]*(state.V[l]*state.V[l]')*channel.H[k,j]')
            Base.LinAlg.BLAS.herk!(Phi_perfect.uplo, 'N', complex(1.), channel.H[k,j]*state.V[l], complex(1.), Phi_perfect.S)
            l != k && Base.LinAlg.BLAS.herk!(Psi_imperfect.uplo, 'N', complex(1.), channel.H[k,j]*state.V[l], complex(1.), Psi_imperfect.S)
        end; end
        for j in setdiff(1:channel.I, coordinators); for l in served_MS_ids(j, assignment)
            Base.LinAlg.BLAS.herk!(Phi_perfect.uplo, 'N', complex(1.), channel.H[k,j]*state.V[l], complex(1.), Phi_perfect.S)
            l != k && robustness && (Psi_imperfect += Hermitian(complex(channel.large_scale_fading_factor[k,j]^2*Ps[j]*eye(channel.Ns[k]))))
        end; end

        # Zero-forcing receiver
        Psi_imperfect_eigen = eigfact(Psi_imperfect, 1:ds[k])
        state.U[k] = Psi_imperfect_eigen.vectors

        # True MMSE receiver and MMSE weight
        F = channel.H[k,i]*state.V[k]
        Ummse = Phi_perfect\F
        state.W[k] = Hermitian((eye(ds[k]) - Ummse'*F)\eye(ds[k]))
    end; end
end

function update_BSs!(state::IntraclusterLeakageMinimizationState,
    channel::SinglecarrierChannel, Ps, assignment, aux_params, robustness)

    ds = [ size(state.W[k], 1) for k = 1:channel.K ]

    for i in active_BSs(assignment)
        coordinators = coordinated_MS_ids(i, assignment)

        # Covariance
        Gamma = Hermitian(complex(zeros(channel.Ms[i], channel.Ms[i])))
        for j = 1:channel.I; for l in served_MS_ids(j, assignment)
            if l in coordinators
                #Gamma += Hermitian(channel.H[k,i]'*(state.U[k]*state.U[k]')*channel.H[k,i])
                Base.LinAlg.BLAS.herk!(Gamma.uplo, 'N', complex(1.), channel.H[l,i]'*state.U[l], complex(1.), Gamma.S)
            else
                if robustness
                    # Using vecnorm(.)^2 may generate NaNs here!
                    Gamma += Hermitian(complex(channel.large_scale_fading_factor[l,i]^2*trace(state.U[l]'*state.U[l])*eye(channel.Ms[i])))
                    # Breaking version: Gamma += Hermitian(complex(channel.large_scale_fading_factor[l,i]^2*vecnorm(state.U[l])^2*eye(channel.Ms[i])))
                end
            end
        end; end

        # Precoders for all served users
        served = served_MS_ids(i, assignment)
        Nserved = length(served)
        for k in served
            Delta = Hermitian(
                Base.LinAlg.BLAS.herk!(Gamma.uplo, 'N', complex(-1.), channel.H[k,i]'*state.U[k], complex(1.), copy(Gamma.S)),
                Gamma.uplo)

            # Precoder
            Delta_eigen = eigfact(Delta, 1:ds[k])
            state.V[k] = sqrt(Ps[i]/(Nserved*ds[k]))*Delta_eigen.vectors
        end
    end
end
