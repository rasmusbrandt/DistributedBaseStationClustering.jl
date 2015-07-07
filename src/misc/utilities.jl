##########################################################################
# Longterm utilities (throughputs, rates and prelogs) for clustering.

# This function calculates the longterm throughputs of the users, given a
# certain network partition (cluster). The rates are the spectral efficiencies
# achievable when IA feasibility is accounted for, and can be -Inf, 0, or
# non-negative depending on the settings. Finally, the throughputs are the
# spectral efficiences when both IA feasibility and the pre-log factor due to
# overhead is accounted for.
function longterm_throughputs(channel, network, partition)
    rates_cluster_sdma, rates_network_sdma = longterm_rates(channel, network, partition)
    prelogs_cluster_sdma, prelogs_network_sdma = longterm_prelogs(network, partition)

    throughputs_cluster_sdma = Diagonal(prelogs_cluster_sdma)*rates_cluster_sdma
    throughputs_network_sdma = Diagonal(prelogs_network_sdma)*rates_network_sdma

    throughputs = throughputs_cluster_sdma .+ throughputs_network_sdma

    return (throughputs, (throughputs_cluster_sdma, throughputs_network_sdma), (rates_cluster_sdma, rates_network_sdma), (prelogs_cluster_sdma, prelogs_network_sdma))
end

# Returns the spectral efficiencies when SDMA is perfomed orthogonally over
# clusters and when SDMA is performed over the network.
function longterm_rates(channel, network, partition)
    I = get_no_BSs(network); K = get_no_MSs(network)
    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    ds = get_no_streams(network); max_d = maximum(ds)
    assignment = get_assignment(network)
    aux_params = get_aux_assignment_params(network)
    IA_infeasible_negative_inf_throughput = aux_params["IA_infeasible_negative_inf_throughput"]

    rates_cluster_sdma = zeros(Float64, K, max_d)
    rates_network_sdma = zeros(Float64, K, max_d)

    for block in partition.blocks
        if is_IA_feasible(network, block)
            # Find out-of-cluster interferers
            intercluster_interferers = setdiff(IntSet(1:I), block.elements) # setdiff is efficient if both arguments are IntSets.
            for i in block.elements
                served = served_MS_ids(i, assignment); Nserved = length(served)
                for k in served
                    desired_power = channel.large_scale_fading_factor[k,i]*channel.large_scale_fading_factor[k,i]*(Ps[i]/(Nserved*ds[k])) # don't user ^2 for performance reasons

                    # Rates when SDMA is perfomed orthogonally over clusters
                    rho = desired_power/sigma2s[k]
                    rates_cluster_sdma[k,1:ds[k]] = longterm_rate(rho)

                    # Rates when SDMA is performed over entire network
                    int_noise_power = sigma2s[k]
                    for j in intercluster_interferers
                        int_noise_power += channel.large_scale_fading_factor[k,j]*channel.large_scale_fading_factor[k,j]*Ps[j] # don't user ^2 for performance reasons
                    end
                    rho = desired_power/int_noise_power
                    rates_network_sdma[k,1:ds[k]] = longterm_rate(rho)
                end
            end
        else
            for i in block.elements; for k in served_MS_ids(i, assignment)
                IA_infeasible_negative_inf_throughput ? rates_network_sdma[k,1:ds[k]] = -Inf : rates_network_sdma[k,1:ds[k]] = 0
            end; end
        end
    end

    return rates_cluster_sdma, rates_network_sdma
end

# Returns prelog factors for the two cases.
function longterm_prelogs(network, partition)
    I = get_no_BSs(network); K = get_no_MSs(network)
    Ns = get_no_MS_antennas(network); Ms = get_no_BS_antennas(network)
    ds = get_no_streams(network)
    assignment = get_assignment(network)
    aux_params = get_aux_assignment_params(network)

    Lcoh = get_aux_network_param(network, "num_coherence_symbols")
    alpha_network_sdma = get_aux_network_param(network, "alpha_network_sdma")
    if !(0 <= alpha_network_sdma <= 1)
        Lumberjack.error("Incorrect alpha_network_sdma. Must be between 0 and 1.")
    end
    alpha_cluster_sdma = 1 - alpha_network_sdma

    # number of symbols owned per BS in the cluster_sdma regime (fair split)
    # we don't floor this to an integer, because we only use this value in the
    # prelog calculation.
    Lorth = alpha_cluster_sdma*Lcoh
    Lorth_per_BS = Lorth/I

    prelogs_cluster_sdma = zeros(Float64, K)
    prelogs_network_sdma = zeros(Float64, K)

    for block in partition.blocks
        cluster_size = length(block.elements)
        Lcsi = CSI_acquisition_symbol_overhead(block, Ns, Ms, ds, assignment) # Symbols needed for CSI acquisition in this cluster
        cluster_sdma_prelog = alpha_cluster_sdma*(cluster_size/I)*max(0., 1 - Lcsi/(cluster_size*Lorth_per_BS))
        for i in block.elements; for k in served_MS_ids(i, assignment)
            prelogs_cluster_sdma[k] = cluster_sdma_prelog

            # Cannot perform network SDMA if the CSI acquisition is infeasible
            if cluster_sdma_prelog == 0.
                prelogs_network_sdma[k] = 0.
            else
                prelogs_network_sdma[k] = alpha_network_sdma
            end
        end; end
    end

    return prelogs_cluster_sdma, prelogs_network_sdma
end

# Calculates longterm rate (or bound thereof) as function of rho
function longterm_rate(rho; bound::Symbol=:none)
    if bound == :upper
        return log2(1 + rho)
    elseif bound == :lower
        return 0.5*log2(1 + 2*rho)
    else
        rho_r = 1/rho
        return (1/log(2))*exp(rho_r + log(expint(rho_r))) # exp(rho_r)*E1(rho_r) can become Inf*0 = NaN
    end
end

# Calculates the number of symbol intervals needed for CSI acquisition in one cluster.
function CSI_acquisition_symbol_overhead(block, Ns, Ms, ds, assignment)
    # First term(DL channel training)
    sum_M = 0
    for i in block.elements
        sum_M += Ms[i]
    end

    # Other terms
    sum_N = 0; quad_sum_M = 0; sum_d = 0
    for i in block.elements
        for k in served_MS_ids(i, assignment)
            # UL channel training
            sum_N += Ns[k]

            # Each MS feeds back all its channels. The feedback matrix thus
            # requires sum_M symbol intervals, assuming that M >= N (see notes).
            quad_sum_M += sum_M

            # DL effective channel training.
            sum_d += ds[k]
        end
    end
    
    return sum_M + sum_N + quad_sum_M + sum_d
end
