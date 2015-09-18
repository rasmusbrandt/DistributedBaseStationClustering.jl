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
    U::Array{Matrix{Complex128},1}
    V::Array{Matrix{Complex128},1}

    # See IntraclusterWMMSEState for description of these.
    E_full::Array{Hermitian{Complex128},1}
    E_partial::Array{Hermitian{Complex128},1}
    E_LB::Array{Hermitian{Complex128},1}
end

NaiveIntraclusterLeakageMinimization(channel, network) =
    IntraclusterLeakageMinimization(channel, network, robustness=false)
RobustIntraclusterLeakageMinimization(channel, network) =
    IntraclusterLeakageMinimization(channel, network, robustness=true)

function IntraclusterLeakageMinimization(channel, network; robustness::Bool=true)
    assignment = get_assignment(network)

    K = get_num_MSs(network)
    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    ds = get_num_streams(network); max_d = maximum(ds)
    alphas = get_user_priorities(network)
    aux_params = get_aux_precoding_params(network)

    clustering_type = get_aux_assignment_param(network, "clustering_type")

    state = IntraclusterLeakageMinimizationState(
        Array(Matrix{Complex128}, K),
        initial_precoders(channel, Ps, sigma2s, ds, assignment, aux_params),
        Array(Hermitian{Complex128}, K),
        Array(Hermitian{Complex128}, K),
        Array(Hermitian{Complex128}, K)
    )
    objective = Float64[]
    weighted_logdet_rates_full = Array(Float64, K, max_d, aux_params["max_iters"])
    weighted_logdet_rates_partial = Array(Float64, K, max_d, aux_params["max_iters"])
    weighted_logdet_rates_LB = Array(Float64, K, max_d, aux_params["max_iters"])
    allocated_power = Array(Float64, K, max_d, aux_params["max_iters"])

    iters = 0; conv_crit = Inf
    while iters < aux_params["max_iters"]
        update_MSs!(state, channel, Ps, sigma2s, assignment,
            clustering_type, robustness)
        iters += 1

        # Results after this iteration
        weighted_logdet_rates_full[:,:,iters], weighted_logdet_rates_partial[:,:,iters], weighted_logdet_rates_LB[:,:,iters] = calculate_weighted_logdet_rates(state, alphas)
        push!(objective, sum(weighted_logdet_rates_full[:,:,iters]))
        allocated_power[:,:,iters] = calculate_allocated_power(state)

        # Check convergence
        if iters >= 2
            conv_crit = abs(objective[end] - objective[end-1])/abs(objective[end-1])
            if conv_crit < aux_params["stop_crit"]
                Lumberjack.debug("IntraclusterLeakageMinimization converged.",
                    @compat Dict(
                        :num_iters => iters,
                        :final_objective => objective[end],
                        :conv_crit => conv_crit,
                        :stop_crit => aux_params["stop_crit"],
                        :max_iters => aux_params["max_iters"]))
                break
            end
        end

        # Begin next iteration, unless the loop will end
        if iters < aux_params["max_iters"]
            update_BSs!(state, channel, Ps, assignment, aux_params,
                clustering_type, robustness)
        end
    end
    if iters == aux_params["max_iters"]
        Lumberjack.debug("IntraclusterLeakageMinimization did NOT converge.",
            @compat Dict(
                :num_iters => iters,
                :final_objective => objective[end],
                :conv_crit => conv_crit,
                :stop_crit => aux_params["stop_crit"],
                :max_iters => aux_params["max_iters"]))
    end

    results = PrecodingResults()
    if aux_params["output_protocol"] == :all_iterations
        results["objective"] = objective
        results["weighted_logdet_rates_full"] = weighted_logdet_rates_full
        results["weighted_logdet_rates_partial"] = weighted_logdet_rates_partial
        results["weighted_logdet_rates_LB"] = weighted_logdet_rates_LB
        results["allocated_power"] = allocated_power
    elseif aux_params["output_protocol"] == :final_iteration
        results["objective"] = objective[iters]
        results["weighted_logdet_rates_full"] = weighted_logdet_rates_full[:,:,iters]
        results["weighted_logdet_rates_partial"] = weighted_logdet_rates_partial[:,:,iters]
        results["weighted_logdet_rates_LB"] = weighted_logdet_rates_LB[:,:,iters]
        results["allocated_power"] = allocated_power[:,:,iters]
    end
    return results
end

function update_MSs!(state::IntraclusterLeakageMinimizationState,
    channel::SinglecarrierChannel, Ps, sigma2s, assignment,
    clustering_type, robustness)

    ds = [ size(state.V[k], 2) for k = 1:channel.K ]

    for i = 1:channel.I; for k in served_MS_ids(i, assignment)
        coordinators = coordinated_BS_ids(k, assignment)

        # Covariances
        Phi_full           = Hermitian(complex(sigma2s[k]*eye(channel.Ns[k]))) # both intercluster and intracluster channels
        Phi_partial_naive  = Hermitian(complex(sigma2s[k]*eye(channel.Ns[k]))) # intercluster instantaneous CSI-R; intracluster ignored
        Phi_partial_robust = Hermitian(complex(sigma2s[k]*eye(channel.Ns[k]))) # intercluster instantaneous CSI-R; intracluster statistical CSI-R
        Psi_full           = Hermitian(complex(sigma2s[k]*eye(channel.Ns[k]))) # both intercluster and intracluster channels
        Psi_partial_naive  = Hermitian(complex(sigma2s[k]*eye(channel.Ns[k]))) # intercluster instantaneous CSI-R; intracluster ignored
        Psi_partial_robust = Hermitian(complex(sigma2s[k]*eye(channel.Ns[k]))) # intercluster instantaneous CSI-R; intracluster statistical CSI-R

        # Intercluster CSI-R
        for j in coordinators; for l in served_MS_ids(j, assignment)
            Fkjl = channel.H[k,j]*state.V[l]
            Base.LinAlg.BLAS.herk!(Phi_full.uplo, 'N', complex(1.), Fkjl, complex(1.), Phi_full.S)
            Base.LinAlg.BLAS.herk!(Phi_partial_naive.uplo, 'N', complex(1.), Fkjl, complex(1.), Phi_partial_naive.S)
            Base.LinAlg.BLAS.herk!(Phi_partial_robust.uplo, 'N', complex(1.), Fkjl, complex(1.), Phi_partial_robust.S)

            if l != k
                Base.LinAlg.BLAS.herk!(Psi_full.uplo, 'N', complex(1.), Fkjl, complex(1.), Psi_full.S)
                Base.LinAlg.BLAS.herk!(Psi_partial_naive.uplo, 'N', complex(1.), Fkjl, complex(1.), Psi_partial_naive.S)
                Base.LinAlg.BLAS.herk!(Psi_partial_robust.uplo, 'N', complex(1.), Fkjl, complex(1.), Psi_partial_robust.S)
            end
        end; end

        # Intracluster CSI-R (only applicable if we are using spectrum_sharing!)
        if clustering_type == :spectrum_sharing
            for j in setdiff(IntSet(1:channel.I), coordinators); for l in served_MS_ids(j, assignment)
                Base.LinAlg.BLAS.herk!(Phi_full.uplo, 'N', complex(1.), channel.H[k,j]*state.V[l], complex(1.), Phi_full.S)
                Phi_partial_robust += Hermitian(complex(channel.large_scale_fading_factor[k,j]^2*Ps[j]*eye(channel.Ns[k])))

                if l != k
                    Base.LinAlg.BLAS.herk!(Psi_full.uplo, 'N', complex(1.), channel.H[k,j]*state.V[l], complex(1.), Psi_full.S)
                    Psi_partial_robust += Hermitian(complex(channel.large_scale_fading_factor[k,j]^2*Ps[j]*eye(channel.Ns[k])))
                end
            end; end
        end

        # Desired channel
        Fiki = channel.H[k,i]*state.V[k]

        # Utility-optimal receiver
        if robustness
            Psi_partial_robust_eigen = eigfact(Psi_partial_robust, 1:ds[k])
            state.U[k] = Psi_partial_robust_eigen.vectors
        else
            Psi_partial_naive_eigen = eigfact(Psi_partial_naive, 1:ds[k])
            state.U[k] = Psi_partial_naive_eigen.vectors
        end

        # "Robust" equation solving for potentially singular effective covariance matrix
        robust_solve(A, B) = try; A\B; catch e; (if isa(e, Base.LinAlg.SingularException); pinv(A)*B; end); end

        # Full CSI used for receive filtering. Intracluster CSI tracked.
        # (This is an achievable rate.)
        U_full = Phi_full\Fiki
        G_full = U_full'*Fiki # effective channel after receive filtering
        Kappa_full = U_full'*Phi_full*U_full # (true) covariance after receive filtering
        S_full = robust_solve(Kappa_full, G_full) # MMSE filter after "original" receive filtering
        state.E_full[k] = Hermitian(eye(ds[k]) - G_full'*S_full)

        # Partial CSI used for receive filtering. Intracluster CSI tracked.
        # (This is an achievable rate.)
        U_partial = state.U[k]
        G_partial = U_partial'*Fiki # effective channel after receive filtering
        Kappa_partial = U_partial'*Phi_full*U_partial # (true) covariance after receive filtering
        S_partial = robust_solve(Kappa_partial, G_partial) # MMSE filter after "original" receive filtering
        state.E_partial[k] = Hermitian(eye(ds[k]) - G_partial'*S_partial)

        # Partial CSI used for receive filtering. Intracluster CSI _not_ tracked.
        # (This is a rate bound)
        U_bound = state.U[k]
        G_bound = U_bound'*Fiki # effective channel after receive filtering
        Kappa_bound = U_bound'*Phi_partial_robust*U_bound # (bound) covariance after receive filtering
        S_bound = robust_solve(Kappa_bound, G_bound) # MMSE filter after "original" receive filtering
        state.E_LB[k] = Hermitian(eye(ds[k]) - G_bound'*S_bound)
    end; end
end

function update_BSs!(state::IntraclusterLeakageMinimizationState,
    channel::SinglecarrierChannel, Ps, assignment, aux_params,
    clustering_type, robustness)

    ds = [ size(state.V[k], 2) for k = 1:channel.K ]

    for i in active_BSs(assignment)
        coordinators = coordinated_MS_ids(i, assignment)

        # Covariance
        Gamma = Hermitian(complex(zeros(channel.Ms[i], channel.Ms[i])))
        for j = 1:channel.I; for l in served_MS_ids(j, assignment)
            if l in coordinators
                #Gamma += Hermitian(channel.H[k,i]'*(state.U[k]*state.U[k]')*channel.H[k,i])
                Base.LinAlg.BLAS.herk!(Gamma.uplo, 'N', complex(1.), channel.H[l,i]'*state.U[l], complex(1.), Gamma.S)
            else
                # Take out-of-cluster interference into account if we are using
                # spectrum sharing, and we actually want to be robust against
                # this type of generated interference.
                if clustering_type == :spectrum_sharing && robustness
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
