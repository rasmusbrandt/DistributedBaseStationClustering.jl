##########################################################################
# Optimal base station clustering based on exhaustive search over the
# utilities proposed in the paper
#
# Chen, Cheng, "Clustering for Interference Alignment in Multiuser
# Interference Network," IEEE Trans. Vehicular Technology, vol. 63, no. 6,
# pp. 2613-2624, July 2014, doi: 10.1109/TVT.2013.2292897

function Chen2014_LinearObj_ExhaustiveSearch(channel, network)
    I = get_no_BSs(network); K = get_no_MSs(network)
    aux_params = get_aux_assignment_params(network)

    # Consistency check
    if I != K
        Lumberjack.error("Chen2014_LinearObj_ExhaustiveSearch can only handle I = K scenarios.")
    end
    if I >= 12
        Lumberjack.warn("Chen2014_LinearObj_ExhaustiveSearch will be slow since I = $I.")
    end

    # Perform cell selection
    LargeScaleFadingCellAssignment!(channel, network)

    # Get W matrix
    W = Chen2014_W_matrix(channel, network)

    # Exhaustive search over partitions
    best_objective = 0.
    best_partition = Partition([0:(I-1)])
    num_sum_throughput_calculations = 0
    for partition in PartitionIterator(I)
        num_sum_throughput_calculations += 1

        # Check that IA is feasible for this cluster structure. Note that this
        # means that Chen2014_LinearObj_ExhaustiveSearch cannot handle situations where
        # IA infeasible blocks are turned off, e.g. when the aux_assignment_param
        # IA_infeasible_negative_inf_throughput is set to false.
        if is_IA_feasible(network, partition)
            # Calculate objective
            objective = 0.
            for block in partition.blocks
                for i in block.elements; for j in block.elements
                    objective += W[i,j]
                end; end
            end

            if objective >= best_objective
                best_objective = objective
                best_partition = partition
            end
        end
    end
    throughputs, _, _, prelogs = longterm_throughputs(channel, network, best_partition)
    a = restricted_growth_string(best_partition)
    Lumberjack.info("Chen2014_LinearObj_ExhaustiveSearch finished.",
        { :sum_throughput => sum(throughputs),
          :a => a }
    )

    # Store prelogs for precoding
    set_aux_network_param!(network, prelogs[1], "prelogs_cluster_sdma")
    set_aux_network_param!(network, prelogs[2], "prelogs_network_sdma")

    # Store cluster assignment together with existing cell assignment
    temp_cell_assignment = get_assignment(network)
    network.assignment = Assignment(network.assignment.cell_assignment, cluster_assignment_matrix(network, best_partition))

    # Return results
    results = AssignmentResults()
    results["throughputs"] = throughputs
    results["a"] = a
    results["num_clusters"] = 1 + maximum(a)
    results["avg_cluster_size"] = avg_cluster_size(a)
    results["num_sum_throughput_calculations"] = num_sum_throughput_calculations
    results["Chen2014_objective"] = best_objective
    return results
end

function Chen2014_kmeans(channel, network)
    I = get_no_BSs(network); K = get_no_MSs(network)
    aux_params = get_aux_assignment_params(network)

    # Consistency check
    if I != K
        Lumberjack.error("Chen2014_kMeans can only handle I = K scenarios.")
    end

    # Get M and N
    require_equal_no_BS_antennas(network)
    require_equal_no_MS_antennas(network)
    M = get_no_BS_antennas(network)[1]
    N = get_no_MS_antennas(network)[1]
    Lmax = M + N - 1

    # Perform cell selection
    LargeScaleFadingCellAssignment!(channel, network)

    # Cluster assignment matrix
    partition_matrix = eye(Int, I, I)

    # Get W matrix
    W = Chen2014_W_matrix(channel, network)

    # k-means clustering on W
    N_A = iceil(I/Lmax)
    if N_A == 1
        # Grand coalition feasible, so no need to perform kmeans
        partition_matrix[:] = 1
    else
        W_eigen = eigfact(W)
        largest_eigenvalues = sortperm(W_eigen.values, rev=true)
        X = W_eigen.vectors[:,largest_eigenvalues[1:N_A]]
        kmeans_clusters = Clustering.kmeans(X', N_A)

        # Convert to partition matrix
        for c = 1:N_A
            ks = find(kmeans_clusters.assignments .== c)
            partition_matrix[ks,ks] = 1
        end
    end

    # Get final throughputs
    partition = Partition(partition_matrix)
    throughputs, _, _, prelogs = longterm_throughputs(channel, network, partition)
    a = restricted_growth_string(partition)
    Lumberjack.info("Chen2014_kMeans finished.",
        { :sum_throughput => sum(throughputs),
          :a => a }
    )

    # Store prelogs for precoding
    set_aux_network_param!(network, prelogs[1], "prelogs_cluster_sdma")
    set_aux_network_param!(network, prelogs[2], "prelogs_network_sdma")

    # Store cluster assignment together with existing cell assignment
    temp_cell_assignment = get_assignment(network)
    network.assignment = Assignment(network.assignment.cell_assignment, cluster_assignment_matrix(network, partition))

    # Return results
    results = AssignmentResults()
    results["throughputs"] = throughputs
    results["a"] = a
    results["num_clusters"] = 1 + maximum(a)
    results["avg_cluster_size"] = avg_cluster_size(a)
    return results
end

function Chen2014_W_matrix(channel, network)
    I = get_no_BSs(network)
    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    assignment = get_assignment(network)

    W = zeros(Float64, I, I)
    for i = 1:I; for j = 1:I
        if i == j
            continue
        end

        k = collect(served_MS_ids(i, assignment))[1]
        desired_power = channel.large_scale_fading_factor[k,i]^2*Ps[i]*channel.large_scale_fading_factor[k,j]^2*Ps[j]
        int_noise_power = sigma2s[k] + channel.large_scale_fading_factor[k,j]^2*Ps[j] + channel.large_scale_fading_factor[k,i]^2*Ps[i]
        W[i,j] = log2(1 + desired_power/int_noise_power)

        l = collect(served_MS_ids(j, assignment))[1]
        desired_power = channel.large_scale_fading_factor[l,j]^2*Ps[j]*channel.large_scale_fading_factor[l,i]^2*Ps[i]
        int_noise_power = sigma2s[k] + channel.large_scale_fading_factor[l,i]^2*Ps[i] + channel.large_scale_fading_factor[l,j]^2*Ps[j]
        W[i,j] += log2(1 + desired_power/int_noise_power)
    end; end

    return W
end
