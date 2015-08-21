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
    I = get_num_BSs(network); K = get_num_MSs(network)
    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    ds = get_num_streams(network); max_d = maximum(ds)
    assignment = get_assignment(network)

    rates_cluster_sdma = zeros(Float64, K, max_d); rates_network_sdma = zeros(Float64, K, max_d)
    for block in partition.blocks
        # If the cluster is not IA feasible, then neither cluster SDMA nor
        # network SDMA will be feasible, and the rates are zero.
        # Note that we _do not_ turn off the BSs that are in IA infeasible
        # clusters; they are thus radiating interference, but no usefuls signal.
        if is_IA_feasible(network, block)
            intercluster_interferers = setdiff(IntSet(1:I), block.elements) # setdiff is efficient if both arguments are IntSets.
            for i in block.elements
                served = served_MS_ids(i, assignment); Nserved = length(served)
                for k in served
                    desired_power = channel.large_scale_fading_factor[k,i]*channel.large_scale_fading_factor[k,i]*(Ps[i]/(Nserved*ds[k])) # don't user ^2 for performance reasons

                    # Rates when SDMA is perfomed orthogonally over clusters
                    rho = desired_power/sigma2s[k]
                    rates_cluster_sdma[k,1:ds[k]] = exp_times_E1(rho)

                    # Rates when SDMA is performed over entire network
                    int_noise_power = sigma2s[k]
                    for j in intercluster_interferers
                        int_noise_power += channel.large_scale_fading_factor[k,j]*channel.large_scale_fading_factor[k,j]*Ps[j] # don't user ^2 for performance reasons
                    end
                    rho = desired_power/int_noise_power
                    rates_network_sdma[k,1:ds[k]] = exp_times_E1(rho)
                end
            end
        end
    end

    return rates_cluster_sdma, rates_network_sdma
end

# Returns prelog factors for the two cases.
function longterm_prelogs(network, partition)
    I = get_num_BSs(network); K = get_num_MSs(network)
    Ns = get_num_MS_antennas(network); Ms = get_num_BS_antennas(network)
    ds = get_num_streams(network)
    assignment = get_assignment(network)

    num_coherence_symbols = get_aux_network_param(network, "num_coherence_symbols")
    beta_network_sdma = get_aux_network_param(network, "beta_network_sdma")
    if !(0 <= beta_network_sdma <= 1)
        Lumberjack.error("Incorrect beta_network_sdma. Must be between 0 and 1.")
    end
    beta_cluster_sdma = 1 - beta_network_sdma

    # Number of symbols for the cluster SDMA phase. We don't floor this
    # to an integer, because we only use this value in the prelog calculation.
    num_symbols_cluster_sdma = beta_cluster_sdma*num_coherence_symbols

    prelogs_cluster_sdma = zeros(Float64, K); prelogs_network_sdma = zeros(Float64, K)
    for block in partition.blocks
        # Calculate cluster SDMA prelog based on CSI acquisition feedback overhead.
        cluster_size = length(block.elements)
        num_CSI_acquisition_symbols = CSI_acquisition_symbol_overhead(block, Ns, Ms, ds, assignment)
        cluster_sdma_prelog = beta_cluster_sdma*max(0., cluster_size/I - num_CSI_acquisition_symbols/num_symbols_cluster_sdma)
        for i in block.elements; for k in served_MS_ids(i, assignment)
            prelogs_cluster_sdma[k] = cluster_sdma_prelog

            # If the CSI acquisition is unbearably high, we cannot
            # do network SDMA for this cluster.
            if cluster_sdma_prelog == 0.
                prelogs_network_sdma[k] = 0.
            else
                prelogs_network_sdma[k] = beta_network_sdma
            end
        end; end
    end

    return prelogs_cluster_sdma, prelogs_network_sdma
end

# Calculates longterm rate (or bound thereof) as function of rho
function exp_times_E1(rho; bound::Symbol=:none)
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
