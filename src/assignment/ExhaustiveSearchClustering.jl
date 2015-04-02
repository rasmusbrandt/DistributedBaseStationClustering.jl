# Exhaustive search over all set partitions with the
# utility model from utilities.jl
function ExhaustiveSearchClustering(channel, network)
    I = get_no_BSs(network); K = get_no_MSs(network)
    d_max = maximum(get_no_streams(network))
    if I > 12
        Lumberjack.warn("ExhaustiveSearchClustering will be slow since I = $I.")
    end

    # Perform cell selection
    LargeScaleFadingCellAssignment!(channel, network)

    # Exhaustive search over all partitions
    no_iters = 0; no_utility_calculations = 0
    best_objective = 0.; best_utilities = Array(Float64, K, d_max)
    best_partition = Partition()
    for partition in PartitionIterator(I)
        no_iters += 1

        # Calculate utilities
        utilities, _ = longterm_utilities(channel, network, partition)
        no_utility_calculations += 1

        objective = sum(utilities)
        if objective > best_objective
            best_objective = objective
            best_utilities = utilities
            best_partition = partition
        end
    end
    a = restricted_growth_string(best_partition)
    Lumberjack.info("ExhaustiveSearchClustering finished.",
        { :sum_utility => best_objective,
          :a => a,
          :no_iters => no_iters,
          :no_utility_calculations => no_utility_calculations }
    )

    # Store cluster assignment together with existing cell assignment
    temp_cell_assignment = get_assignment(network)
    network.assignment = Assignment(temp_cell_assignment.cell_assignment, cluster_assignment_matrix(network, best_partition))

    # Return results
    results = AssignmentResults()
    results["utilities"] = best_utilities
    results["a"] = a
    results["no_iters"] = no_iters
    results["no_utility_calculations"] = no_utility_calculations
    return results
end
