##########################################################################
# Non-cooperative clustering, i.e. all coalitions are singleton.

function NoClustering(channel, network)
    I = get_num_BSs(network); K = get_num_MSs(network)
    aux_params = get_aux_assignment_params(network)

    # Perform cell selection
    LargeScaleFadingCellAssignment!(channel, network)
    temp_assignment = get_assignment(network)

    a = collect(0:(I-1))
    partition = Partition(a)
    utilities, alphas, _ = longterm_utilities(channel, network, partition)
    Lumberjack.info("NoClustering finished.",
        @compat Dict(
            :sum_utility => sum(utilities),
            :a => a))

    # Store alphas as user priorities for precoding, if desired
    if aux_params["apply_overhead_prelog"]
        set_user_priorities!(network, alphas)
    end

    # Store cluster assignment together with existing cell assignment
    network.assignment = Assignment(temp_assignment.cell_assignment, cluster_assignment_matrix(network, partition))

    # Return results
    results = AssignmentResults()
    results["utilities"] = utilities
    results["a"] = a
    results["alphas"] = alphas
    results["num_clusters"] = 1 + maximum(a)
    results["num_utility_calculations"] = K
    results["num_longterm_rate_calculations"] = K
    return results
end
