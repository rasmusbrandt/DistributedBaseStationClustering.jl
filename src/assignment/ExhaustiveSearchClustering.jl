function ExhaustiveSearchClustering(channel, network)
    K = get_no_MSs(network)
    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    ds = get_no_streams(network)

    partitions = all_partitions(1:K)

    # Exhaustive search over partitions
    rates = Array(Float64, K)
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
            intercluster_interferers = setdiff(1:K, block.elements)
            for k in block.elements
                desired_power = channel.large_scale_fading_factor[k,k]^2*Ps[k]
                int_noise_power = sigma2s[k]
                for l in intercluster_interferers
                    int_noise_power += channel.large_scale_fading_factor[k,l]^2*Ps[l]
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
    assignment_matrix = partition_to_assignment_matrix(best_partition, K)
    # println(assignment_matrix)

    assign_cells_by_id!(network)
    network.assignment = Assignment(network.assignment.cell_assignment, assignment_matrix)
end
