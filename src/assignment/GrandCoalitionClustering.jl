# Helper function that gives the grand coalition, i.e. that all BSs
# cooperate within one huge cluster.
function GrandCoalitionClustering(channel, network)
    I = get_no_BSs(network); K = get_no_MSs(network)

    # Perform cell selection
    LargeScaleFadingCellAssignment!(channel, network)

    a = zeros(Int, I)
    utilities, _ = longterm_utilities(channel, network, Partition(a))
    Lumberjack.info("GrandCoalitionClustering finished.", { :sum_utility => sum(utilities), :a => a })

    # Store cluster assignment together with existing cell assignment
    temp_assignment = get_assignment(network)
    network.assignment = Assignment(temp_assignment.cell_assignment, ones(K, I))

    # Return results
    results = AssignmentResults()
    results["utilities"] = utilities
    results["a"] = a
    results["no_iters"] = 1
    results["no_utility_calculations"] = 1
    return results
end
