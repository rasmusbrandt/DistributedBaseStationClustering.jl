##########################################################################
# Utilities for clustering and coalitional formation

# This function calculates the longterm utilities of the users, given a
# certain network partition (cluster). The utopian_rates are the spectral
# efficiences achievable disregarding IA feasibility and overheads, and are
# thus upper bounds. The rates are the spectral efficiencies achievable
# when IA feasibility is accounted for, and can be -Inf, 0, or non-negative
# depending on the settings. Finally, the throughputs are the spectral
# efficiences when both IA feasibility and the pre-log factor due to
# overhead is accounted for. The alphas are the pre-log factors.
#
# If the spectral efficiency when IA is not feasible is 0, this means that
# the corresponding cluster is turned off. This does however _not_ mean that
# the overhead pre-log factor is changed. We still count the overhead from
# CSI acquisition, and in the case of orthogonal clustering the temporal
# degrees of freedom provided by each BS in the clusters, when calculating
# the overheads. This could/should be generalized to the case when some
# clusters are turned off, and thus giving up their temporal degrees of
# freedom and CSI acquisition overhead to the other clusters. In the current
# implementation, we are not taking this into consideration however.
#
# Notice: this function only returns a bound of the rates for now,
# since I don't if the E1 exponential integral has been implemented
# by anyone in Julia yet.
function longterm_utilities(channel, network, partition)
    I = get_no_BSs(network); K = get_no_MSs(network)
    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    ds = get_no_streams(network); max_d = maximum(ds)
    assignment = get_assignment(network)
    aux_params = get_aux_assignment_params(network)

    @defaultize_param! aux_params "IntraclusterWMMSE:bisection_Gamma_cond" 1e10

    utopian_rates = zeros(Float64, K, max_d) # raw spectral efficiency upper bound, disregarding IA feasiblility and model applicability
    rates = zeros(Float64, K, max_d) # raw spectral efficiency, zero or -Inf if IA not feasible
    throughputs = zeros(Float64, K, max_d) # spectral efficiency incl. overhead pre-log factor
    alphas = ones(Float64, K) # pre-log factor due to overhead

    # Both rate and overhead calculations depend on what type of
    # clustering that is used.
    if aux_params["clustering_type"] == :orthogonal
        # In orthogonal clustering, only the direct channel affects the rate.
        for block in partition.blocks
            # But the pre-log factors can be different between coalitions.
            alpha = orthogonal_prelog_factor(network, block)

            # Check IA feasibility for this block
            IA_feas = is_IA_feasible(network, block)

            # Calculate rates
            for i in block.elements
                served = served_MS_ids(i, assignment); Nserved = length(served)
                for k in served
                    # Retain overhead pre-log factors, to be used in the precoding
                    alphas[k] = alpha

                    # Rates without interference, assuming IA feasibility
                    desired_power = channel.large_scale_fading_factor[k,i]^2*(Ps[i]/(Nserved*ds[k]))
                    rho = desired_power/sigma2s[k]

                    utopian_rates[k,1:ds[k]] = 0.5log(1 + 2rho) # This is a lower bound

                    if IA_feas
                        # These rates are achievable using IA
                        rates[k,:] = utopian_rates[k,:]
                        throughputs[k,:] = alphas[k]*utopian_rates[k,:]
                    else
                        # Not feasible for IA. Set rates for this block
                        # as 0 or -Inf, depending on given parameter.
                        if aux_params["IA_infeasible_utility_inf"]
                            rates[k,1:ds[k]] = -Inf; throughputs[k,1:ds[k]] = -Inf
                        else
                            rates[k,:] = 0; throughputs[k,:] = 0
                        end
                    end
                end
            end
        end
    elseif aux_params["clustering_type"] == :spectrum_sharing
        # Only the BSs that are in the partition will be active in generating
        # interference to the other 
        active_BSs = IntSet()
        for block in partition.blocks
            union!(active_BSs, block.elements)
        end

        # Calculate rates for all MSs in clusters
        for block in partition.blocks
            # Check IA feasibility for this block
            IA_feas = is_IA_feasible(network, block)

            # Find out-of-cluster interferers
            intercluster_interferers = setdiff(active_BSs, block.elements) # setdiff is efficient if both arguments are IntSets.
            for i in block.elements
                served = served_MS_ids(i, assignment); Nserved = length(served)
                for k in served
                    # Rates without interference, assuming IA feasibility
                    desired_power = channel.large_scale_fading_factor[k,i]*channel.large_scale_fading_factor[k,i]*(Ps[i]/(Nserved*ds[k])) # don't user ^2 for performance reasons
                    int_noise_power = sigma2s[k]
                    for j in intercluster_interferers
                        int_noise_power += channel.large_scale_fading_factor[k,j]*channel.large_scale_fading_factor[k,j]*Ps[j] # don't user ^2 for performance reasons
                    end
                    rho = desired_power/int_noise_power

                    utopian_rates[k,1:ds[k]] = 0.5log(1 + 2rho) # This is a lower bound

                    if IA_feas
                        # These rates are achievable using IA
                        rates[k,:] = utopian_rates[k,:]
                        throughputs[k,:] = utopian_rates[k,:] # will be multiplied by alpha later
                    else
                        # Not feasible for IA. Set rates for this block
                        # as 0 or -Inf, depending on given parameter.
                        if aux_params["IA_infeasible_utility_inf"]
                            rates[k,1:ds[k]] = -Inf; throughputs[k,1:ds[k]] = -Inf
                        else
                            rates[k,:] = 0; throughputs[k,:] = 0
                        end
                    end
                end
            end
        end

        # For spectrum sharing, the pre-log factor is identical over coalitions.
        alpha = spectrum_sharing_prelog_factor(network, partition)
        throughputs = alpha*rates
        alphas *= alpha
    end

    # By having apply_overhead_prelog as an assignment parameter, we can easily
    # run the simulations with and without overhead prelog applied.
    if aux_params["apply_overhead_prelog"]
        return throughputs, alphas, utopian_rates
    else
        return rates, alphas, utopian_rates
    end
end

# Pre-log factor for spectrum sharing clustering. All BSs already share the
# entire coherence time, so adding BSs to coalitions monotonically increases
# the overhead.
function spectrum_sharing_prelog_factor(network, partition)
    Tc = get_aux_network_param(network, "no_coherence_symbols")

    # Sum of Lps (implemented as loop, since it seems to be the fastest way)
    Lp_sum = 0
    for block in partition.blocks
        Lp_sum += CSI_acquisition_symbol_overhead(network, block)
    end

    return max(0., 1 - Lp_sum/Tc)
end

# Pre-log factor for orthogonal clustering. Each BS brings its share (1/I) of
# the coherence time to the coalition it joins, but also contributes extra
# overhead due to the CSI acquisition.
function orthogonal_prelog_factor(network, block)
    I = get_no_BSs(network)
    Tc = get_aux_network_param(network, "no_coherence_symbols")
    return max(0., length(block)/I - CSI_acquisition_symbol_overhead(network, block)/Tc)
end

# Calculates the number of symbol intervals needed for CSI acquisition. This
# quantity is denoted with $L_p^\text{CSI}$ in the notes.
function CSI_acquisition_symbol_overhead(network, block)
    Ns = get_no_MS_antennas(network); Ms = get_no_BS_antennas(network)
    ds = get_no_streams(network)
    assignment = get_assignment(network)

    # First term in Lp (DL channel training)
    sum_M = 0
    for i in block.elements
        sum_M += Ms[i]
    end

    # Other terms in Lp
    sum_N = 0; quad_sum_M = 0; sum_d = 0
    for i in block.elements
        for k in served_MS_ids(i, assignment)
            # UL channel training
            sum_N += Ns[k]

            # Each MS feeds back all its channels. The feedback matrix thus
            # requires sum_M symbol intervals, assuming that M > N (see notes).
            quad_sum_M += sum_M

            # DL effective channel training.
            sum_d += ds[k]
        end
    end
    
    return sum_M + sum_N + quad_sum_M + sum_d
end
