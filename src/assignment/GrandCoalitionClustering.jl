function GrandCoalitionClustering(channel, network)
    I = get_no_BSs(network); K = get_no_MSs(network)

    # Perform cell selection
    LargeScaleFadingCellAssignment!(channel, network)

    a = zeros(Int, I)
    Lumberjack.info("GrandCoalitionClustering finished.", { :sum_rate => sum(longterm_throughputs(channel, network, Partition(a))), :a => a })

    # Store cluster assignment together with existing cell assignment
    temp_assignment = get_assignment(network)
    network.assignment = Assignment(temp_assignment.cell_assignment, ones(K, I))
end
