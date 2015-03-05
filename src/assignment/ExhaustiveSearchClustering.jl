function ExhaustiveSearchClustering(channel, network)
    I = get_no_BSs(network); K = get_no_MSs(network)
    ds = get_no_streams(network)

    # Perform cell selection
    LargeScaleFadingCellAssignment!(channel, network)
    temp_cell_assignment = get_assignment(network)

    # Exhaustive search over all partitions
    rates = Array(Float64, K)
    best_partition = Partition(); best_sum_rate = 0.
    for partition in PartitionIterator(I)
        # Check that IA is feasible for this cluster structure
        if is_IA_feasible(partition, channel.Ns, channel.Ms, ds, temp_cell_assignment)
            # Calculate rates
            rates = longterm_cluster_IA_rates(channel, network, partition, temp_cell_assignment)

            sum_rate = sum(rates)
            if sum_rate > best_sum_rate
                best_sum_rate = sum_rate
                best_partition = partition
            end
        end
    end

    # Build cluster assignment matrix
    cluster_assignment_matrix = partition_to_cluster_assignment_matrix(best_partition, K, I, temp_cell_assignment)

    # Store cluster assignment together with existing cell assignment
    network.assignment = Assignment(temp_cell_assignment.cell_assignment, cluster_assignment_matrix)
end
