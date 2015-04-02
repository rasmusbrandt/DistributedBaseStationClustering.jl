# A simple greedy clustering method based on path loss and transmit powers.
# In each step, the method finds the strongest interfering link, and matches
# the corresponding cells into a cluster, if that is IA feasible. Thus, the
# bound on the cluster sizes will be given by IA feasibility, and not by the
# overhead pre-log factor.
function GreedyClustering(channel, network)
    I = get_no_BSs(network); K = get_no_MSs(network)
    aux_params = get_aux_assignment_params(network)

    # Consistency check
    if aux_params["IA_infeasible_negative_inf_utility"] == false
        Lumberjack.info("GreedyClustering only finds solutions where all clusters are IA feasible. IA_infeasible_negative_inf_utility is set to false, which means that the other methods might find solutions where some blocks are turned off due to IA infeasibility.")
    end

    # Perform cell selection
    LargeScaleFadingCellAssignment!(channel, network)
    temp_cell_assignment = get_assignment(network)

    # Cluster assignment matrix
    partition_matrix = eye(Int, I, I)

    # Clustering metric
    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    F = zeros(Float64, I, I)
    for i = 1:I; for j = 1:I
        if i == j
            F[i,j] = -Inf
            continue
        end

        # Sum interference from BS j to MSs served by BS i
        F[i,j] = sum([ log2((channel.large_scale_fading_factor[k,j]^2)*Ps[j]/sigma2s[k]) for k in served_MS_ids(i, temp_cell_assignment) ])
    end; end

    # Greedily build clusters based on strongest sum interference between cells
    no_iters = 0
    while !all(F .== -Inf)
        no_iters += 1

        # Find strongest interfering link that is still active
        _, idx = findmax(F)
        i, j = ind2sub((I, I), idx)
        i_cluster = find(partition_matrix[i,:] .== 1)

        # Assign to cluster
        partition_matrix[i_cluster,j] = 1; partition_matrix[j,i_cluster] = 1

        # Check that IA is feasible for this cluster structure. Note that this
        # means that GreedyClustering cannot handle situations where
        # IA infeasible blocks are turned off, e.g. when the aux_assignment_param
        # IA_infeasible_negative_inf_utility is set to false.
        if is_IA_feasible(network, Partition(partition_matrix))
            # Fix BS j to this cluster
            F[:,j] = -Inf

            if length(i_cluster) == 1
                # We have effectively put BS i in this cluster as well,
                # without it really being aware. Do not try to put BS i
                # in another cluster.
                F[:,i] = -Inf
            end
        else
            # This was not a feasible cluster, undo the assignment.
            partition_matrix[i_cluster,j] = 0; partition_matrix[j,i_cluster] = 0

            # Do not try to add BS j to the cluster of BS i again
            F[i_cluster,j] = -Inf
        end
    end
    partition = Partition(partition_matrix)
    utilities, _ = longterm_utilities(channel, network, partition)
    objective = sum(utilities)
    Lumberjack.info("GreedyClustering finished.",
        { :sum_utility => objective,
          :a => restricted_growth_string(partition_matrix),
          :no_iters => no_iters }
    )

    # Store cluster assignment together with existing cell assignment
    network.assignment = Assignment(temp_cell_assignment.cell_assignment, cluster_assignment_matrix(network, partition))

    # Return results
    results = AssignmentResults()
    results["utilities"] = utilities
    results["a"] = restricted_growth_string(partition_matrix)
    return results
end
