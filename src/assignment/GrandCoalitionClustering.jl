function GrandCoalitionClustering(channel, network)
    I = get_no_BSs(network)

    # Perform cell selection
    LargeScaleFadingCellAssignment!(channel, network)
    temp_assignment = get_assignment(network)

    # Store cluster assignment together with existing cell assignment
    network.assignment = Assignment(temp_assignment.cell_assignment, ones(I, I))
end
