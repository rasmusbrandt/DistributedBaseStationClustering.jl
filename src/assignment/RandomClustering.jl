function RandomClustering(channel, network)
    I = get_no_BSs(network)
    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    ds = get_no_streams(network)

    # Perform cell selection
    LargeScaleFadingCellAssignment!(channel, network)
    temp_assignment = get_assignment(network)

    partitions = all_partitions(1:I)

    # Pick random partition
    selected_partition = partitions[int(floor(I*rand()) + 1)]

    # Build cluster assignment matrix
    assignment_matrix = partition_to_assignment_matrix(selected_partition, I)

    # Store cluster assignment together with existing cell assignment
    network.assignment = Assignment(temp_assignment.cell_assignment, assignment_matrix)
end
