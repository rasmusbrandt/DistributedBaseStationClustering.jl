##########################################################################
# Optimal base station clustering based on exhaustive search.
#
# All possible restricted growth strings (and thus set partitions) are
# enumerated, and the best (in the utilities.jl utilities sense) is picked.

function ExhaustiveSearchClustering(channel, network)
    I = get_num_BSs(network); K = get_num_MSs(network)
    d_max = maximum(get_num_streams(network))
    aux_params = get_aux_assignment_params(network)

    # Warn if this will be slow...
    if I >= 12
        Lumberjack.warn("ExhaustiveSearchClustering will be slow since I = $I.")
    end

    # Perform cell selection
    LargeScaleFadingCellAssignment!(channel, network)

    # Exhaustive search over all partitions
    num_utility_calculations = 0
    num_longterm_rate_calculations = sum([ binomial(I,i) for i = 1:I ]) # number of possible distinct clusters whose members need to calculate their longterm rates
    best_objective = 0.; best_utilities = Array(Float64, K, d_max)
    best_alphas = Array(Float64, K); best_partition = Partition(collect(0:(I-1)))
    for partition in PartitionIterator(I)
        num_utility_calculations += K

        # Calculate utilities
        utilities, alphas, _ = longterm_utilities(channel, network, partition)

        objective = sum(utilities)
        if objective >= best_objective
            best_objective = objective
            best_utilities = utilities
            best_alphas = alphas
            best_partition = partition
        end
    end
    a = restricted_growth_string(best_partition)
    Lumberjack.info("ExhaustiveSearchClustering finished.",
        @compat Dict(
            :sum_utility => best_objective,
            :a => a))

    # Store alphas as user priorities for precoding, if desired
    if aux_params["apply_overhead_prelog"]
        set_user_priorities!(network, best_alphas)
    end

    # Store cluster assignment together with existing cell assignment
    temp_cell_assignment = get_assignment(network)
    network.assignment = Assignment(temp_cell_assignment.cell_assignment, cluster_assignment_matrix(network, best_partition))

    # Return results
    results = AssignmentResults()
    results["utilities"] = best_utilities
    results["a"] = a
    results["alphas"] = best_alphas
    results["num_clusters"] = 1 + maximum(a)
    results["num_utility_calculations"] = num_utility_calculations
    results["num_longterm_rate_calculations"] = num_longterm_rate_calculations
    return results
end
