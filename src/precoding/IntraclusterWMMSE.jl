immutable IntraclusterWMMSEState
    U::Array{Matrix{Complex128},1}
    W::Array{Hermitian{Complex128},1} # for rate calculations
    Z::Array{Hermitian{Complex128},1} # MSE weights actually used
    V::Array{Matrix{Complex128},1}
end

NaiveIntraclusterWMMSE(channel, network) = IntraclusterWMMSE(channel, network, robustness=false)
RobustIntraclusterWMMSE(channel, network) = IntraclusterWMMSE(channel, network, robustness=true)

function IntraclusterWMMSE(channel, network; robustness=true)
    assignment = get_assignment(network)

    K = get_no_MSs(network)
    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    ds = get_no_streams(network)

    aux_params = get_aux_precoding_params(network)
    @defaultize_param! aux_params "IntraclusterWMMSE:bisection_Gamma_cond" 1e10
    @defaultize_param! aux_params "IntraclusterWMMSE:bisection_singular_Gamma_mu_lower_bound" 1e-14
    @defaultize_param! aux_params "IntraclusterWMMSE:bisection_max_iters" 5e1
    @defaultize_param! aux_params "IntraclusterWMMSE:bisection_tolerance" 1e-3

    state = IntraclusterWMMSEState(
        Array(Matrix{Complex128}, K),
        unity_MSE_weights(ds),
        unity_MSE_weights(ds),
        initial_precoders(channel, Ps, sigma2s, ds, assignment, aux_params))
    objective = Float64[]
    utilities = Array(Float64, K, maximum(ds), aux_params["max_iters"])
    logdet_rates = Array(Float64, K, maximum(ds), aux_params["max_iters"])
    MMSE_rates = Array(Float64, K, maximum(ds), aux_params["max_iters"])
    allocated_power = Array(Float64, K, maximum(ds), aux_params["max_iters"])

    iters = 0; conv_crit = Inf
    while iters < aux_params["max_iters"]
        update_MSs!(state, channel, Ps, sigma2s, assignment, robustness)
        iters += 1

        # Results after this iteration
        utilities[:,:,iters] = calculate_utilities(state)
        push!(objective, sum(utilities[:,:,iters]))
        logdet_rates[:,:,iters] = calculate_logdet_rates(state)
        MMSE_rates[:,:,iters] = calculate_MMSE_rates(state)
        allocated_power[:,:,iters] = calculate_allocated_power(state)

        # Check convergence
        if iters >= 2
            conv_crit = abs(objective[end] - objective[end-1])/abs(objective[end-1])
            if conv_crit < aux_params["stop_crit"]
                Lumberjack.debug("IntraclusterWMMSE converged.",
                    { :no_iters => iters, :final_objective => objective[end],
                      :conv_crit => conv_crit, :stop_crit => aux_params["stop_crit"],
                      :max_iters => aux_params["max_iters"] })
                break
            end
        end

        # Begin next iteration, unless the loop will end
        if iters < aux_params["max_iters"]
            update_BSs!(state, channel, Ps, assignment, aux_params, robustness)
        end
    end
    if iters == aux_params["max_iters"]
        Lumberjack.debug("IntraclusterWMMSE did NOT converge.",
            { :no_iters => iters, :final_objective => objective[end],
              :conv_crit => conv_crit, :stop_crit => aux_params["stop_crit"],
              :max_iters => aux_params["max_iters"] })
    end

    results = PrecodingResults()
    if aux_params["output_protocol"] == :all_iterations
        results["objective"] = objective
        results["utilities"] = utilities
        results["logdet_rates"] = logdet_rates
        results["MMSE_rates"] = MMSE_rates
        results["allocated_power"] = allocated_power
    elseif aux_params["output_protocol"] == :final_iteration
        results["objective"] = objective[iters]
        results["utilities"] = utilities[:,:,iters]
        results["logdet_rates"] = logdet_rates[:,:,iters]
        results["MMSE_rates"] = MMSE_rates[:,:,iters]
        results["allocated_power"] = allocated_power[:,:,iters]
    end
    return results
end

function update_MSs!(state::IntraclusterWMMSEState,
    channel::SinglecarrierChannel, Ps, sigma2s, assignment, robustness)

    ds = [ size(state.Z[k], 1) for k = 1:channel.K ]

    for i = 1:channel.I; for k in served_MS_ids(i, assignment)
        coordinators = coordinated_BS_ids(k, assignment)

        # Covariances
        Phi_imperfect = Hermitian(complex(sigma2s[k]*eye(channel.Ns[k])))
        Phi_perfect = Hermitian(complex(sigma2s[k]*eye(channel.Ns[k])))
        for j = 1:channel.I
            if j in coordinators
                for l in served_MS_ids(j, assignment)
                    #Phi += Hermitian(channel.H[k,j]*(state.V[l]*state.V[l]')*channel.H[k,j]')
                    Base.LinAlg.BLAS.herk!(Phi_imperfect.uplo, 'N', complex(1.), channel.H[k,j]*state.V[l], complex(1.), Phi_imperfect.S)
                    Base.LinAlg.BLAS.herk!(Phi_perfect.uplo, 'N', complex(1.), channel.H[k,j]*state.V[l], complex(1.), Phi_perfect.S)
                end
            else
                for l in served_MS_ids(j, assignment)
                    if robustness
                        Phi_imperfect += Hermitian(complex(channel.large_scale_fading_factor[k,j]^2*Ps[j]*eye(channel.Ns[k])))
                    end
                    Base.LinAlg.BLAS.herk!(Phi_perfect.uplo, 'N', complex(1.), channel.H[k,j]*state.V[l], complex(1.), Phi_perfect.S)
                end
            end
        end

        # Intracluster receiver and MSE weight
        F = channel.H[k,i]*state.V[k]
        state.U[k] = Phi_imperfect\F
        state.Z[k] = Hermitian((eye(ds[k]) - state.U[k]'*F)\eye(ds[k]))

        # True MMSE receiver and MMSE weight
        Ummse = Phi_perfect\F
        state.W[k] = Hermitian((eye(ds[k]) - Ummse'*F)\eye(ds[k]))
    end; end
end

function update_BSs!(state::IntraclusterWMMSEState,
    channel::SinglecarrierChannel, Ps, assignment, aux_params, robustness)

    for i in active_BSs(assignment)
        coordinators = coordinated_MS_ids(i, assignment)

        # Virtual uplink covariance
        Gamma = Hermitian(complex(zeros(channel.Ms[i],channel.Ms[i])))
        for j = 1:channel.I; for l in served_MS_ids(j, assignment)
            if l in coordinators
                Gamma += Hermitian(channel.H[l,i]'*(state.U[l]*state.Z[l]*state.U[l]')*channel.H[l,i])
            else
                if robustness
                    Gamma += Hermitian(complex(channel.large_scale_fading_factor[l,i]^2*vecnorm(state.U[l])^2*eye(channel.Ms[i])))
                end
            end
        end; end

        # Find optimal Lagrange multiplier
        mu_star, Gamma_eigen =
            optimal_mu(i, Gamma, state, channel, Ps, assignment, aux_params)

        # Precoders (reuse EVD)
        for k in served_MS_ids(i, assignment)
            state.V[k] = Gamma_eigen.vectors*Diagonal(1./(Gamma_eigen.values .+ mu_star))*Gamma_eigen.vectors'*channel.H[k,i]'*state.U[k]*state.Z[k]
        end
    end
end

function optimal_mu(i, Gamma, state::IntraclusterWMMSEState,
    channel::SinglecarrierChannel, Ps, assignment, aux_params)

    # Build bisector function
    bis_M = Hermitian(complex(zeros(channel.Ms[i], channel.Ms[i])))
    for k in served_MS_ids(i, assignment)
        #bis_M += Hermitian(channel.H[k,i]'*(state.U[k]*(state.Z[k]*state.Z[k])*state.U[k]')*channel.H[k,i])
        Base.LinAlg.BLAS.herk!(bis_M.uplo, 'N', complex(1.), channel.H[k,i]'*state.U[k]*state.Z[k], complex(1.), bis_M.S)
    end
    Gamma_eigen = eigfact(Gamma)
    bis_JMJ_diag = abs(diag(Gamma_eigen.vectors'*bis_M*Gamma_eigen.vectors))
    f(mu) = sum(bis_JMJ_diag./((Gamma_eigen.values .+ mu).*(Gamma_eigen.values .+ mu)))

    # mu lower bound
    if abs(maximum(Gamma_eigen.values))/abs(minimum(Gamma_eigen.values)) < aux_params["IntraclusterWMMSE:bisection_Gamma_cond"]
        # Gamma is invertible
        mu_lower = 0
    else
        mu_lower = aux_params["IntraclusterWMMSE:bisection_singular_Gamma_mu_lower_bound"]
    end

    if f(mu_lower) <= Ps[i]
        # No bisection needed
        return mu_lower, Gamma_eigen
    else
        # mu upper bound
        mu_upper = sqrt(channel.Ms[i]/Ps[i]*maximum(bis_JMJ_diag)) - abs(minimum(Gamma_eigen.values))
        if f(mu_upper) > Ps[i]
            error("Power bisection: infeasible mu upper bound.")
        end

        no_iters = 0
        while no_iters < aux_params["IntraclusterWMMSE:bisection_max_iters"]
            conv_crit = (Ps[i] - f(mu_upper))/Ps[i]

            if conv_crit < aux_params["IntraclusterWMMSE:bisection_tolerance"]
                break
            else
                mu = (1/2)*(mu_lower + mu_upper)

                if f(mu) < Ps[i]
                    # New point feasible, replace upper point
                    mu_upper = mu
                else
                    # New point not feasible, replace lower point
                    mu_lower = mu
                end
            end

            no_iters += 1
        end

        if no_iters == aux_params["IntraclusterWMMSE:bisection_max_iters"]
            println("Power bisection: reached max iterations.")
        end

        # The upper point is always feasible, therefore we use it
        return mu_upper, Gamma_eigen
    end
end

function calculate_utilities(state::IntraclusterWMMSEState)
    K = length(state.Z)
    ds = Int[ size(state.Z[k], 1) for k = 1:K ]; max_d = maximum(ds)

    utilities = zeros(Float64, K, max_d)
    for k = 1:K; if ds[k] > 0
        # W is p.d., so we should only get abs eigenvalues. Numerically we may
        # get some imaginary noise however. Also, numerically the eigenvalues
        # may be less than 1, so we need to handle that to not get negative
        # rates.
        r = log2(max(1, abs(eigvals(state.Z[k]))))

        if ds[k] < max_d
            utilities[k,:] = cat(1, r, zeros(Float64, max_d - ds[k]))
        else
            utilities[k,:] = r
        end
    end; end

    return utilities
end
