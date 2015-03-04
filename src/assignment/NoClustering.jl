function NoClustering(channel, network)
    I = get_no_BSs(network); K = get_no_MSs(network)

    # Perform cell selection
    LargeScaleFadingCellAssignment!(channel, network)
    temp_assignment = get_assignment(network)

    cluster_assignment_matrix = zeros(Int, K, I)
    for i = 1:I; for k in served_MS_ids(i, temp_assignment)
        cluster_assignment_matrix[k,i] = 1
    end; end

    # Store cluster assignment together with existing cell assignment
    network.assignment = Assignment(temp_assignment.cell_assignment, cluster_assignment_matrix)
end
