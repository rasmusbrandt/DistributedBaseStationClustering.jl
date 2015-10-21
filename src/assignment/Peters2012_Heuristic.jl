##########################################################################
# Greedy base station clustering from
#
# Peters, Heath, "User Partitioning for Less Overhead in MIMO
# Interference Channels", IEEE Trans. WC, vol. 11, no. 2, pp. 592-603, 2012.

function Peters2012_Heuristic(channel, network)
    I = get_num_BSs(network); K = get_num_MSs(network)
    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    aux_params = get_aux_assignment_params(network)

    # Ensure that the network is symmetric (needed for DoF calculation)
    require_equal_num_BS_antennas(network); M = get_num_BS_antennas(network)[1]
    require_equal_num_MS_antennas(network); N = get_num_MS_antennas(network)[1]
    require_equal_num_streams(network); d = get_num_streams(network)[1]

    # Perform cell selection
    LargeScaleFadingCellAssignment!(channel, network)
    temp_cell_assignment = get_assignment(network)

    # Find DoF optimal number of partitions P
    DoFo = zeros(Float64, I) # d-tilde in paper
    for L = 1:I # cluster size
        tmp_block = Block(IntSet(1:L)) # pseudo-block for overhead pre-log factor calc.
        DoFo[L] = orthogonal_prelog_factor(network, tmp_block)*floor((N + M)*L/(L+1))
    end
    _, idx = findmax(DoFo)
    P = ceil(Int, I/(1:I)[idx]) # number of partitions

    # Greedily build clusters based on rate approximation
    unclustered_BSs = IntSet(1:I) # K_A in paper
    clusters = [ IntSet() for p = 1:P ] # K_P in paper

    num_iters = 0
    while length(unclustered_BSs) > 0
        num_iters += 1

        # Build clustering metric
        rate_approx = zeros(Float64, I, P) # already clustered BSs will have zero metric
        for i in unclustered_BSs
            served = served_MS_ids(i, temp_cell_assignment)
            Nserved = length(served)
            for p in 1:P
                # pseudo-block for overhead pre-log factor calc.
                tmp_block = Block(union(clusters[p], IntSet(i)))

                # Peters only supports the IC. We make the obvious extension
                # to the IBC here however.
                for k in served
                    # We use the pathloss information rather than the actual channel strength.
                    # This is contrary to Peters metric, which is the Frobenius norm of the
                    # instantaneous channel. This modified metric makes sense in our
                    # clustering framework however.
                    rate_approx[i,p] += orthogonal_prelog_factor(network, tmp_block)*d*log2(1 + (channel.large_scale_fading_factor[k,i]^2)*(Ps[i]/Nserved)/sigma2s[k])
                end
            end
        end
        _, idx = findmax(rate_approx)
        i_star, p_star = ind2sub((I, P), idx)

        # Add BS to cluster
        push!(clusters[p_star], i_star)

        # Don't cluster BS again
        delete!(unclustered_BSs, i_star)
    end

    # Build partition from clusters
    partition = Partition()
    for c in clusters
        if length(c) > 0
            push!(partition.blocks, Block(c))
        end
    end

    # Get result
    utilities, alphas, _ = longterm_utilities(channel, network, partition)
    a = restricted_growth_string(partition)
    objective = sum(utilities)
    Lumberjack.info("Peters2010_Heuristic finished.",
        @compat Dict(
            :sum_utility => objective,
            :a => a))

    # Store alphas as user priorities for precoding, if desired
    if aux_params["apply_overhead_prelog"]
        set_user_priorities!(network, alphas)
    end

    # Store cluster assignment together with existing cell assignment
    network.assignment = Assignment(temp_cell_assignment.cell_assignment, cluster_assignment_matrix(network, partition))

    # Return results
    results = AssignmentResults()
    results["utilities"] = utilities
    results["a"] = a
    results["alphas"] = alphas
    results["num_clusters"] = 1 + maximum(a)
    results["num_iters"] = num_iters
    return results
end
