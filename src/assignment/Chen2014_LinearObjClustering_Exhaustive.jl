function Chen2014_LinearObjClustering_Exhaustive(channel, network)
    I = get_no_BSs(network); K = get_no_MSs(network)
    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    ds = get_no_streams(network)

    # Consistency check
    if I != K
        Lumberjack.error("Chen2014_LinearObjClustering can only handle I = K scenarios.")
    end

    # Perform cell selection
    LargeScaleFadingCellAssignment!(channel, network)
    temp_assignment = get_assignment(network)

    W = zeros(Float64, I, I)
    for i = 1:I; for j = 1:I
        if i == j
            continue
        end

        k = collect(served_MS_ids(i, temp_assignment))[1]
        desired_power = channel.large_scale_fading_factor[k,i]^2*Ps[i]*channel.large_scale_fading_factor[k,j]^2*Ps[j]
        int_noise_power = sigma2s[k] + channel.large_scale_fading_factor[k,j]^2*Ps[j] + channel.large_scale_fading_factor[k,i]^2*Ps[i]
        W[i,j] = log2(1 + desired_power/int_noise_power)

        l = collect(served_MS_ids(j, temp_assignment))[1]
        desired_power = channel.large_scale_fading_factor[l,j]^2*Ps[j]*channel.large_scale_fading_factor[l,i]^2*Ps[i]
        int_noise_power = sigma2s[k] + channel.large_scale_fading_factor[l,i]^2*Ps[i] + channel.large_scale_fading_factor[l,j]^2*Ps[j]
        W[i,j] += log2(1 + desired_power/int_noise_power)
    end; end

    partitions = all_partitions(1:I)

    # Exhaustive search over partitions
    objective = 0.
    best_partition = Partition(); best_objective = 0.
    for partition in partitions
        # Check that IA is feasible for this cluster structure
        if is_IA_feasible(partition, channel.Ns, channel.Ms, ds, temp_assignment)
            # Calculate objective
            for block in partition.blocks
                objective = 0.
                for i in block.elements
                    for j in block.elements
                        objective += W[i,j]
                    end
                end
            end

            if objective > best_objective
                best_objective = objective
                best_partition = partition
            end
        end
    end

    # Build cluster assignment matrix
    assignment_matrix = partition_to_assignment_matrix(best_partition, K, I, temp_assignment)

    # Store cluster assignment together with existing cell assignment
    network.assignment = Assignment(network.assignment.cell_assignment, assignment_matrix)
end
