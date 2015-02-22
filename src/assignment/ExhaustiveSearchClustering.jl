function ExhaustiveSearchClustering(channel, network)
    I = get_no_BSs(network); K = get_no_MSs(network)
    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    ds = get_no_streams(network)

    # Perform cell selection
    LargeScaleFadingCellAssignment!(channel, network)
    temp_assignment = get_assignment(network)

    # Exhaustive search over all partitions
    rates = Array(Float64, K)
    partitions = all_partitions(1:I)
    best_partition = Partition(); best_sum_rate = 0.
    for partition in partitions
        # Check that IA is feasible for this cluster structure
        if is_IA_feasible(partition, channel.Ns, channel.Ms, ds, temp_assignment)
            # Calculate rates
            for block in partition.blocks
                intercluster_interferers = setdiff(1:I, block.elements)
                for i in block.elements
                    for k in served_MS_ids(i, temp_assignment)
                        desired_power = channel.large_scale_fading_factor[k,i]^2*Ps[i]
                        int_noise_power = sigma2s[k]
                        for j in intercluster_interferers
                            int_noise_power += channel.large_scale_fading_factor[k,j]^2*Ps[j]
                        end
                        rho = 1 + desired_power/int_noise_power

                        rates[k] = 0.5log(1 + 2rho)
                    end
                end
            end

            sum_rate = sum(rates)
            if sum_rate > best_sum_rate
                best_sum_rate = sum_rate
                best_partition = partition
            end
        end
    end

    # Build cluster assignment matrix
    assignment_matrix = partition_to_assignment_matrix(best_partition, K, I, temp_assignment)

    # Store cluster assignment together with existing cell assignment
    network.assignment = Assignment(temp_assignment.cell_assignment, assignment_matrix)
end
