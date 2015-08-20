##########################################################################
# Optimal base station clustering based on exhaustive search.
#
# All possible restricted growth strings (and thus set partitions) are
# enumerated, and the best (in the utilities.jl throughputs sense) is picked.

function ExhaustiveSearchClustering(channel, network)
    I = get_no_BSs(network); K = get_no_MSs(network)
    d_max = maximum(get_no_streams(network))
    aux_params = get_aux_assignment_params(network)

    # Warn if this will be slow...
    if I >= 12
        Lumberjack.warn("ExhaustiveSearchClustering will be slow since I = $I.")
    end

    # Perform cell selection
    LargeScaleFadingCellAssignment!(channel, network)

    # Exhaustive search over all partitions
    num_sum_throughput_calculations = 0
    best_sum_throughput = 0.; best_throughputs = Array(Float64, K, d_max); local best_throughputs_split
    best_prelogs = (Array(Float64, K), Array(Float64, K)); best_partition = Partition([0:(I-1)])
    for partition in PartitionIterator(I)
        num_sum_throughput_calculations += 1

        # Calculate throughputs
        throughputs, throughputs_split, _, prelogs = longterm_throughputs(channel, network, partition)

        objective = sum(throughputs)
        if objective >= best_sum_throughput
            best_sum_throughput = objective
            best_throughputs = throughputs
            best_throughputs_split = throughputs_split
            best_prelogs = prelogs
            best_partition = partition
        end
    end
    a = restricted_growth_string(best_partition)
    Lumberjack.info("ExhaustiveSearchClustering finished.",
        { :sum_throughput => best_sum_throughput,
          :num_sum_throughput_calculations => num_sum_throughput_calculations,
          :a => a }
    )

    # Store prelogs for precoding
    set_aux_network_param!(network, best_prelogs[1], "prelogs_cluster_sdma")
    set_aux_network_param!(network, best_prelogs[2], "prelogs_network_sdma")

    # Store cluster assignment together with existing cell assignment
    temp_cell_assignment = get_assignment(network)
    network.assignment = Assignment(temp_cell_assignment.cell_assignment, cluster_assignment_matrix(network, best_partition))

    # Return results
    results = AssignmentResults()
    results["throughputs"] = best_throughputs
    results["throughputs_cluster_sdma"] = best_throughputs_split[1]
    results["throughputs_network_sdma"] = best_throughputs_split[2]
    results["a"] = a
    results["num_clusters"] = 1 + maximum(a)
    results["avg_cluster_size"] = avg_cluster_size(a)
    results["num_sum_throughput_calculations"] = num_sum_throughput_calculations
    return results
end
