##########################################################################
# Non-cooperative clustering, i.e. all coalitions are singleton.

function NoClustering(channel, network)
    I = get_no_BSs(network); K = get_no_MSs(network)
    aux_params = get_aux_assignment_params(network)

    # Perform cell selection
    LargeScaleFadingCellAssignment!(channel, network)
    temp_assignment = get_assignment(network)

    cluster_assignment_matrix = zeros(Int, K, I)
    for i = 1:I; for k in served_MS_ids(i, temp_assignment)
        cluster_assignment_matrix[k,i] = 1
    end; end

    a = [0:(I-1)]
    utilities, alphas, _ = longterm_utilities(channel, network, Partition(a))
    Lumberjack.info("NoClustering finished.",
        { :sum_utility => sum(utilities),
          :a => a }
    )

    # Store alphas as user priorities for precoding, if desired
    if aux_params["apply_overhead_prelog"]
        set_user_priorities!(network, alphas)
    end

    # Store cluster assignment together with existing cell assignment
    network.assignment = Assignment(temp_assignment.cell_assignment, cluster_assignment_matrix)

    # Return results
    results = AssignmentResults()
    results["utilities"] = utilities
    results["a"] = a
    results["alphas"] = alphas
    results["no_clusters"] = 1 + maximum(a)
    results["no_utility_calculations"] = K
    results["no_longterm_rate_calculations"] = K
    return results
end
