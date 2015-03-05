function ExhaustiveSearchClustering(channel, network)
    I = get_no_BSs(network); K = get_no_MSs(network)

    # Perform cell selection
    LargeScaleFadingCellAssignment!(channel, network)

    # Exhaustive search over all partitions
    rates = Array(Float64, K)
    best_partition = Partition(); best_sum_rate = 0.
    for partition in PartitionIterator(I)
        # Check that IA is feasible for this cluster structure
        if is_IA_feasible(network, partition)
            # Calculate rates
            rates = longterm_cluster_IA_rates(channel, network, partition)

            sum_rate = sum(rates)
            if sum_rate > best_sum_rate
                best_sum_rate = sum_rate
                best_partition = partition
            end
        end
    end

    # Build cluster assignment matrix
    cluster_assignment_matrix = partition_to_cluster_assignment_matrix(network, best_partition)

    # Store cluster assignment together with existing cell assignment
    temp_cell_assignment = get_assignment(network)
    network.assignment = Assignment(temp_cell_assignment.cell_assignment, cluster_assignment_matrix)
end
