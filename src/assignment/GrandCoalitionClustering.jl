##########################################################################
# Grand coalition base station clustering.
#
# This means that all BSs cooperate with each other. Note that this might
# mean that the throughputs are either zero or -Inf, since the utility model
# is only applicable when IA is feasible.

function GrandCoalitionClustering(channel, network)
    I = get_num_BSs(network); K = get_num_MSs(network)
    aux_params = get_aux_assignment_params(network)

    # Perform cell selection
    LargeScaleFadingCellAssignment!(channel, network)

    a = zeros(Int, I)
    throughputs, throughputs_split, _, prelogs =
        longterm_throughputs(channel, network, Partition(a))
    Lumberjack.info("GrandCoalitionClustering finished.",
        @compat Dict(:sum_throughput => sum(throughputs), :a => a))

    # Store prelogs for precoding
    set_aux_network_param!(network, prelogs[1], "prelogs_cluster_sdma")
    set_aux_network_param!(network, prelogs[2], "prelogs_network_sdma")

    # Store cluster assignment together with existing cell assignment
    temp_assignment = get_assignment(network)
    network.assignment = Assignment(temp_assignment.cell_assignment, ones(K, I))

    # Return results
    results = AssignmentResults()
    results["throughputs"] = throughputs
    results["throughputs_cluster_sdma"] = throughputs_split[1]
    results["throughputs_network_sdma"] = throughputs_split[2]
    results["a"] = a
    results["num_clusters"] = reshape(1 + maximum(a), 1, 1)
    results["avg_cluster_size"] = reshape(avg_cluster_size(a), 1, 1)
    results["num_sum_throughput_calculations"] = reshape(1, 1, 1)
    return results
end
