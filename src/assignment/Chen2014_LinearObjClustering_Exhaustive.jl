function Chen2014_LinearObjClustering_Exhaustive(channel, network)
    K = get_no_MSs(network)
    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    ds = get_no_streams(network)

    W = zeros(Float64, K, K)
    for k = 1:K; for l = 1:K
        if k == l
            continue
        end

        desired_power = channel.large_scale_fading_factor[k,k]^2*Ps[k]*channel.large_scale_fading_factor[k,l]^2*Ps[l]
        int_noise_power = sigma2s[k] + channel.large_scale_fading_factor[k,l]^2*Ps[l] + channel.large_scale_fading_factor[k,k]^2*Ps[k]
        W[k,l] = log2(1 + desired_power/int_noise_power)

        desired_power = channel.large_scale_fading_factor[l,l]^2*Ps[l]*channel.large_scale_fading_factor[l,k]^2*Ps[k]
        int_noise_power = sigma2s[k] + channel.large_scale_fading_factor[l,k]^2*Ps[k] + channel.large_scale_fading_factor[l,l]^2*Ps[l]
        W[k,l] += log2(1 + desired_power/int_noise_power)
    end; end

    partitions = all_partitions(1:K)

    # Exhaustive search over partitions
    objective = 0.
    best_partition = Partition(); best_objective = 0.
    for partition in partitions
        feasible = true

        for block in partition.blocks
            # IA feasibility check for this cluster
            if !feasibility_IA(partition, channel.Ns, channel.Ms, ds)
                feasible = false
                break
            end

            # Calculate objective for this cluster
            objective = 0.
            for k in block.elements
                for l in block.elements
                    objective += W[k,l]
                end
            end
        end

        if feasible
            # println(objective)
            if objective > best_objective
                best_objective = objective
                best_partition = partition
            end
        end
    end

    # Build cluster assignment matrix
    assignment_matrix = partition_to_assignment_matrix(best_partition, K)
    # println(assignment_matrix)

    assign_cells_by_id!(network)
    network.assignment = Assignment(network.assignment.cell_assignment, assignment_matrix)
end
