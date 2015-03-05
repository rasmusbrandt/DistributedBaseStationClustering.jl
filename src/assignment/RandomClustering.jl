function RandomClustering(channel, network)
    I = get_no_BSs(network); K = get_no_MSs(network)
    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    ds = get_no_streams(network)

    # Perform cell selection
    LargeScaleFadingCellAssignment!(channel, network)
    temp_cell_assignment = get_assignment(network)

    # Generate random restricted growth string. See Algorithm H in TAoCP 7.2.1.5.
    a = Array(Int, I)
    b = Array(Int, I)
    for i = 1:I
        if i == 1
            b[1] = 1
            a[1] = 0
        else
            b[i] = 1 + maximum(a[1:(i-1)])
            a[i] = ifloor((b[i] + 1)*rand())
        end
    end

    # Get corresponding partition
    partition = Partition(a)

    # Build cluster assignment matrix
    assignment_matrix = partition_to_cluster_assignment_matrix(partition, K, I, temp_cell_assignment)

    # Store cluster assignment together with existing cell assignment
    network.assignment = Assignment(temp_cell_assignment.cell_assignment, assignment_matrix)
end
