function ExhaustiveSearchClustering(channel, network)
    I = get_no_BSs(network)
    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    ds = get_no_streams(network)

    # Perform cell selection
    assign_cells_by_large_scale_fading!(channel, network)
    temp_assignment = get_assignment(network)

    # Exhaustive search over all partitions
    rates = Array(Float64, I)
    partitions = all_partitions(1:I)
    best_partition = Partition(); best_sum_rate = 0.
    for partition in partitions
        feasible = true

        for block in partition.blocks
            # IA feasibility check for this cluster
            if !feasibility_IA(partition, channel.Ns, channel.Ms, ds)
                feasible = false
                break
            end

            # Calculate rates
            intercluster_interferers = setdiff(1:I, block.elements)
            for i in block.elements
                k = collect(served_MS_ids(i, temp_assignment))[1]
                desired_power = channel.large_scale_fading_factor[k,i]^2*Ps[i]
                int_noise_power = sigma2s[k]
                for j in intercluster_interferers
                    int_noise_power += channel.large_scale_fading_factor[k,j]^2*Ps[j]
                end
                rho = 1 + desired_power/int_noise_power

                rates[k] = 0.5log(1 + 2rho)
            end
        end

        if feasible
            sum_rate = sum(rates)
            # println(sum_rate)
            if sum_rate > best_sum_rate
                best_sum_rate = sum_rate
                best_partition = partition
            end
        end
    end

    # Build cluster assignment matrix
    assignment_matrix = partition_to_assignment_matrix(best_partition, I)

    # Store cluster assignment together with existing cell assignment
    network.assignment = Assignment(temp_assignment.cell_assignment, assignment_matrix)
end
