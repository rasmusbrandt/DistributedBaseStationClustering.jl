function GrandCoalitionClustering(channel, network)
    I = get_no_BSs(network); K = get_no_MSs(network)

    # Perform cell selection
    LargeScaleFadingCellAssignment!(channel, network)

    Lumberjack.info("GrandCoalitionClustering finished.", { :a => zeros(Int, I) })

    # Store cluster assignment together with existing cell assignment
    temp_assignment = get_assignment(network)
    network.assignment = Assignment(temp_assignment.cell_assignment, ones(K, I))
end
