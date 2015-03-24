##########################################################################
# Utilities for clustering and coalitional formation.

# This function calculates the longterm throughputs of the users, given a
# certain network partition (cluster). If a block in the partition is not
# IA feasible, the throughputs are -Inf for those users. The return vector
# can thus contain both Float64s, as well as Infs. It is up to the caller
# to select how to handle this.
#
# Notice: this function only returns a bound of the rates for now,
# since I don't if the E1 exponential integral has been implemented
# by anyone in Julia yet.
function longterm_throughputs(channel, network, partition; apply_overhead_prelog::Bool=true)
    I = get_no_BSs(network); K = get_no_MSs(network)
    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    ds = get_no_streams(network)
    assignment = get_assignment(network)
    aux_params = get_aux_assignment_params(network)

    # Rates are the raw spectral efficients. Throughputs are the
    # spectral efficiencies with overhead pre-log factor applied.
    rates = zeros(Float64, K)

    # Both rate and overhead calculations depend on what type of
    # clustering that is used.
    if aux_params["clustering_type"] == :orthogonal
        # In orthogonal clustering, only the direct channel affects the rate.
        for block in partition.blocks
            # But the pre-log factors can be different between coalitions.
            alpha = orthogonal_prelog_factor(network, block)

            # Calculate rates
            for i in block.elements
                served = served_MS_ids(i, assignment); Nserved = length(served)
                for k in served
                    desired_power = channel.large_scale_fading_factor[k,i]^2*(Ps[i]/(Nserved*ds[k]))
                    rho = desired_power/sigma2s[k]

                    rates[k] = ds[k]*0.5log(1 + 2rho) # This is a lower bound

                    # Calculate throughputs
                    if apply_overhead_prelog
                        throughputs[k] = alpha*rates[k]
                    else
                        throughputs[k] = rates[k]
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
            intercluster_interferers = setdiff(active_BSs, block.elements)
            for i in block.elements
                served = served_MS_ids(i, assignment); Nserved = length(served)
                for k in served
                    desired_power = channel.large_scale_fading_factor[k,i]^2*(Ps[i]/(Nserved*ds[k]))
                    int_noise_power = sigma2s[k] + sum([ channel.large_scale_fading_factor[k,j]^2*Ps[j] for j in intercluster_interferers ])
                    rho = desired_power/int_noise_power

                    rates[k] = ds[k]*0.5log(1 + 2rho) # This is a lower bound
                end
            end
        end

        # For spectrum sharing, the pre-log factor is identical over coalitions.
        if apply_overhead_prelog
            throughputs = spectrum_sharing_prelog_factor(network, partition)*rates
        else
            throughputs = rates
        end
    else
        Lumberjack.error("Incorrect clustering given in auxiliary assignment parameters.")
    end

    return throughputs
end

# The long-term rates are simply the long-term throughputs, when the overhead
# pre-log factors are 1. This is implemented by simply not multiplying with the
# overhead pre-log factors.
longterm_rates(channel, network, partition) =
    longterm_throughputs(channel, network, partition, apply_overhead_prelog=false)

# Pre-log factor for spectrum sharing clustering. All BSs already share the
# entire coherence time, so adding BSs to coalitions monotonically increases
# the overhead.
function spectrum_sharing_prelog_factor(network, partition)
    Tc = get_aux_assignment_param(network, "no_coherence_symbols")
    return max(0, 1 - sum([ CSI_acquisition_symbol_overhead(network, block) for block in partition.blocks ])/Tc)
end

# Pre-log factor for orthogonal clustering. Each BS brings its share (1/I) of
# the coherence time to the coalition it joins, but also contributes extra
# overhead due to the CSI acquisition.
function orthogonal_prelog_factor(network, block)
    I = get_no_BSs(network)
    Tc = get_aux_assignment_param(network, "no_coherence_symbols")
    return max(0, length(block)/I - CSI_acquisition_symbol_overhead(network, block)/Tc)
end

# Calculates the number of symbol intervals needed for CSI acquisition. This
# quantity is denoted with $L_p^\text{CSI}$ in the notes.
function CSI_acquisition_symbol_overhead(network, block)
    Ns = get_no_MS_antennas(network); Ms = get_no_BS_antennas(network)
    ds = get_no_streams(network)
    assignment = get_assignment(network)

    # First term in Lp (DL channel training)
    sum_M = sum(Ms[collect(block.elements)])

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
