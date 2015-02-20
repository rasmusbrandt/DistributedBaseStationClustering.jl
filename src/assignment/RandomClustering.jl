function RandomClustering(channel, network)
    K = get_no_MSs(network)
    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    ds = get_no_streams(network)

    partitions = all_partitions(1:K)

    # Pick random partition
    selected_partition = partitions[int(floor(K*rand()) + 1)]

    # Build cluster assignment matrix
    assignment_matrix = partition_to_assignment_matrix(selected_partition, K)
    println(assignment_matrix)

    assign_cells_by_id!(network)
    network.assignment = Assignment(network.assignment.cell_assignment, assignment_matrix)
end
