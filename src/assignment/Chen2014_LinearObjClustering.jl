function Chen2014_LinearObjClustering_ExhaustiveSearch(channel, network)
    I = get_no_BSs(network); K = get_no_MSs(network)

    # Consistency check
    if I != K
        Lumberjack.error("Chen2014_LinearObjClustering can only handle I = K scenarios.")
    end

    # Perform cell selection
    LargeScaleFadingCellAssignment!(channel, network)

    # Get W matrix
    W = Chen2014_LinearObjClustering_W(channel, network)

    # Exhaustive search over partitions
    best_partition = Partition(); best_objective = 0.
    for partition in PartitionIterator(I)
        # Check that IA is feasible for this cluster structure
        if is_IA_feasible(network, partition)
            # Calculate objective
            objective = 0.
            for block in partition.blocks
                objective = 0.
                for i in block.elements; for j in block.elements
                    objective += W[i,j]
                end; end
            end

            if objective > best_objective
                best_objective = objective
                best_partition = partition
            end
        end
    end
    Lumberjack.info("Chen2014_LinearObjClustering_ExhaustiveSearch finished.",
        { :best_objective => best_objective,
          :best_partition => best_partition }
    )

    # Build cluster assignment matrix
    cluster_assignment_matrix = partition_to_cluster_assignment_matrix(network, best_partition)

    # Store cluster assignment together with existing cell assignment
    temp_cell_assignment = get_assignment(network)
    network.assignment = Assignment(network.assignment.cell_assignment, cluster_assignment_matrix)
end

function Chen2014_LinearObjClustering_W(channel, network)
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
