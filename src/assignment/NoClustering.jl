##########################################################################
# Non-cooperative clustering, i.e. all coalitions are singleton.

function NoClustering(channel, network)
    I = get_no_BSs(network); K = get_no_MSs(network)
    aux_params = get_aux_assignment_params(network)

    # Perform cell selection
    LargeScaleFadingCellAssignment!(channel, network)
    temp_assignment = get_assignment(network)

    a = [0:(I-1)]
    partition = Partition(a)
    throughputs, _, _, prelogs = longterm_throughputs(channel, network, partition)
    Lumberjack.info("NoClustering finished.",
        { :sum_throughput => sum(throughputs),
          :a => a }
    )

    # Store prelogs for precoding
    set_aux_network_param!(network, prelogs[1], "prelogs_cluster_sdma")
    set_aux_network_param!(network, prelogs[2], "prelogs_network_sdma")

    # Store cluster assignment together with existing cell assignment
    network.assignment = Assignment(temp_assignment.cell_assignment, cluster_assignment_matrix(network, partition))

    # Return results
    results = AssignmentResults()
    results["throughputs"] = throughputs
    results["a"] = a
    results["num_clusters"] = 1 + maximum(a)
    results["avg_cluster_size"] = avg_cluster_size(a)
    results["num_sum_throughput_calculations"] = 1
    return results
end
