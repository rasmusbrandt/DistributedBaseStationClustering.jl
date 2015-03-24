function ExhaustiveSearchClustering(channel, network)
    I = get_no_BSs(network); K = get_no_MSs(network)

    # Perform cell selection
    LargeScaleFadingCellAssignment!(channel, network)

    # Exhaustive search over all partitions
    no_iters = 0
    best_partition = Partition(); best_objective = 0.
    for partition in PartitionIterator(I)
        no_iters += 1

        # Check that IA is feasible for this cluster structure
        if is_IA_feasible(network, partition)
            # Calculate rates
            rates = longterm_throughputs(channel, network, partition)

            objective = sum(rates)
            if objective > best_objective
                best_objective = objective
                best_partition = partition
            end
        end
    end
    Lumberjack.info("ExhaustiveSearchClustering finished.",
        { :sum_rate => best_objective,
          :a => restricted_growth_string(best_partition),
          :no_iters => no_iters }
    )

    # Store cluster assignment together with existing cell assignment
    temp_cell_assignment = get_assignment(network)
    network.assignment = Assignment(temp_cell_assignment.cell_assignment, cluster_assignment_matrix(network, best_partition))
end
