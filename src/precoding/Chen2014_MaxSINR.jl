##########################################################################
# Chen2014_MaxSINR
#
# Out-of-cluster interference robustified version of MaxSINR from
# Chen, Cheng, "Clustering for Interference Alignment in Multiuser
# Interference Network", IEEE Trans. VT, vol. 63, no. 6, pp. 2613-2624, 2014.

immutable Chen2014_MaxSINRState
    U::Array{Matrix{Complex128},1}
    W::Array{Hermitian{Complex128},1} # w/ intracluster CSI available
    Z::Array{Hermitian{Complex128},1} # w/o intracluster CSI available
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

    clustering_type = get_aux_assignment_param(network, "clustering_type")

    state = Chen2014_MaxSINRState(
        initial_receivers(channel, Ps, sigma2s, ds, assignment, aux_params),
        Array(Hermitian{Complex128}, K),
        Array(Hermitian{Complex128}, K),
        initial_precoders(channel, Ps, sigma2s, ds, assignment, aux_params)
    )
    objective = Float64[]
    utilities = Array(Float64, K, max_d, aux_params["max_iters"])
    logdet_rates = Array(Float64, K, max_d, aux_params["max_iters"])
    MMSE_rates = Array(Float64, K, max_d, aux_params["max_iters"])
    weighted_logdet_rates = Array(Float64, K, max_d, aux_params["max_iters"])
    weighted_MMSE_rates = Array(Float64, K, max_d, aux_params["max_iters"])
    allocated_power = Array(Float64, K, max_d, aux_params["max_iters"])

    iters = 0; conv_crit = Inf
    while iters < aux_params["max_iters"]
        update_MSs!(state, channel, Ps, sigma2s, assignment,
            clustering_type, robustness)
        iters += 1

        # Results after this iteration
        utilities[:,:,iters] = calculate_utilities(state, alphas)
        push!(objective, sum(utilities[:,:,iters]))
        logdet_rates[:,:,iters] = calculate_logdet_rates(state)
        MMSE_rates[:,:,iters] = calculate_MMSE_rates(state)
        weighted_logdet_rates[:,:,iters] = calculate_weighted_logdet_rates(state, alphas)
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
            update_BSs!(state, channel, Ps, sigma2s, assignment, aux_params,
                clustering_type, robustness)
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
        results["utilities"] = utilities
        results["logdet_rates"] = logdet_rates
        results["MMSE_rates"] = MMSE_rates
        results["weighted_logdet_rates"] = weighted_logdet_rates
        results["weighted_MMSE_rates"] = weighted_MMSE_rates
        results["allocated_power"] = allocated_power
    elseif aux_params["output_protocol"] == :final_iteration
        results["objective"] = objective[iters]
        results["utilities"] = utilities[:,:,iters]
        results["logdet_rates"] = logdet_rates[:,:,iters]
        results["MMSE_rates"] = MMSE_rates[:,:,iters]
        results["weighted_logdet_rates"] = weighted_logdet_rates[:,:,iters]
        results["weighted_MMSE_rates"] = weighted_MMSE_rates[:,:,iters]
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
        Phi_perfect = Hermitian(complex(sigma2s[k]*eye(channel.Ns[k])))
        Phi_imperfect = Hermitian(complex(sigma2s[k]*eye(channel.Ns[k])))
        for j in coordinators; for l in served_MS_ids(j, assignment)
            #Phi += Hermitian(channel.H[k,j]*(state.V[l]*state.V[l]')*channel.H[k,j]')
            Base.LinAlg.BLAS.herk!(Phi_perfect.uplo, 'N', complex(1.), channel.H[k,j]*state.V[l], complex(1.), Phi_perfect.S)
            Base.LinAlg.BLAS.herk!(Phi_imperfect.uplo, 'N', complex(1.), channel.H[k,j]*state.V[l], complex(1.), Phi_imperfect.S)
        end; end

        # Take out-of-cluster interference into account if we are using
        # spectrum sharing clustering. If we are using orthogonal clustering,
        # there is no out-of-cluster interference.
        if clustering_type == :spectrum_sharing
            for j in setdiff(1:channel.I, coordinators); for l in served_MS_ids(j, assignment)
                Base.LinAlg.BLAS.herk!(Phi_perfect.uplo, 'N', complex(1.), channel.H[k,j]*state.V[l], complex(1.), Phi_perfect.S)
                robustness && (Phi_imperfect += Hermitian(complex(channel.large_scale_fading_factor[k,j]^2*Ps[j]*eye(channel.Ns[k]))))
            end; end
        end

        # Per-stream receivers
        for n = 1:ds[k]
            Psi_plus_n = Hermitian(
                Base.LinAlg.BLAS.herk!(Phi_imperfect.uplo, 'N', complex(-1.), channel.H[k,i]*state.V[k][:,n], complex(1.), copy(Phi_imperfect.S)),
                Phi_imperfect.uplo)
            u = Psi_plus_n\channel.H[k,i]*state.V[k][:,n]
            state.U[k][:,n] = u/norm(u,2)
        end

        # Intracluster receiver and MSE weight
        F = channel.H[k,i]*state.V[k]
        U = Phi_imperfect\F
        state.Z[k] = Hermitian((eye(ds[k]) - U'*F)\eye(ds[k]))

        # True MMSE receiver and MMSE weight
        Ummse = Phi_perfect\F
        state.W[k] = Hermitian((eye(ds[k]) - Ummse'*F)\eye(ds[k]))
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

# Weighted rates when the spatial characteristics of the out-of-cluster
# interference is unknown to the receivers, but all intracluster interference
# is perfectly known.
function calculate_utilities(state::Chen2014_MaxSINRState, alphas)
    K = length(state.V)
    ds = Int[ size(state.V[k], 2) for k = 1:K ]; max_d = maximum(ds)

    utilities = zeros(Float64, K, max_d)
    for k = 1:K; if ds[k] > 0
        # W is p.d., so we should only get abs eigenvalues. Numerically we may
        # get some imaginary noise however. Also, numerically the eigenvalues
        # may be less than 1, so we need to handle that to not get negative
        # rates.
        r = alphas[k]*log2(max(1, abs(eigvals(state.Z[k]))))

        if ds[k] < max_d
            utilities[k,:] = cat(1, r, zeros(Float64, max_d - ds[k]))
        else
            utilities[k,:] = r
        end
    end; end

    return utilities
end
