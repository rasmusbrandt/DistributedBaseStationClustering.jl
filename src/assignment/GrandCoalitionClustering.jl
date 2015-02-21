function GrandCoalitionClustering(channel, network)
    I = get_no_BSs(network)

    # Perform cell selection
    assign_cells_by_large_scale_fading!(channel, network)
    temp_assignment = get_assignment(network)

    # Store cluster assignment together with existing cell assignment
    network.assignment = Assignment(temp_assignment.cell_assignment, ones(I, I))
end
