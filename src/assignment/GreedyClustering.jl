##########################################################################
# Greedy base station clustering.
#
# A simple greedy clustering method based on path loss and transmit powers.
# In each step, the method finds the strongest interfering link, and matches
# the corresponding cells into a cluster, if that is IA feasible. This is
# done either by putting the offending BS into the cluster of the victim
# BS (GreedyClustering_Single), or by trying to merge the respective
# clusters, if possible (GreedyClustering_Multiple). The bound on the
# cluster sizes will be given by IA feasibility, and not by the overhead
# pre-log factor. If IA_infeasible_negative_inf_utility is set to false,
# it will mean that other methods might find solutions where some BSs
# are put in clusters that are turned off due to IA infeasibility. These
# methods will not be able to do that.

GreedyClustering_Single(channel, network) =
    GreedyClustering(channel, network, merge_multiple=false)
GreedyClustering_Multiple(channel, network) =
    GreedyClustering(channel, network, merge_multiple=true)

function GreedyClustering(channel, network; merge_multiple::Bool=false)
    I = get_num_BSs(network); K = get_num_MSs(network)
    aux_params = get_aux_assignment_params(network)

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
    num_utility_calculations = 0
    num_longterm_rate_calculations = 0
    while !all(F .== -Inf)
        # Find strongest interfering link that is still active
        _, idx = findmax(F)
        i, j = ind2sub((I, I), idx)

        if merge_multiple
            # Join clusters of BS i and BS j
            i_cluster = find(partition_matrix[i,:] .== 1)
            j_cluster = find(partition_matrix[j,:] .== 1)
            new_partition_matrix = copy(partition_matrix)
            new_partition_matrix[i_cluster,j_cluster] = 1
            new_partition_matrix[j_cluster,i_cluster] = 1

            # Check IA feasibility for this new cluster. (Note that this
            # means that GreedyClustering cannot handle situations where
            # IA infeasible blocks are turned off, e.g. when the aux_assignment_param
            # IA_infeasible_negative_inf_utility is set to false.)
            num_utility_calculations += K
            num_longterm_rate_calculations += length(i_cluster) + length(j_cluster)
            if is_IA_feasible(network, Partition(new_partition_matrix))
                partition_matrix = new_partition_matrix
            end

             # Never consider this link again
            F[i,j] = -Inf
        else
            # Assign BS j to cluster of BS i
            i_cluster = find(partition_matrix[i,:] .== 1)
            new_partition_matrix = copy(partition_matrix)
            new_partition_matrix[i_cluster,j] = 1; new_partition_matrix[j,i_cluster] = 1

            # Check IA feasibility for this new cluster.
            num_utility_calculations += K
            num_longterm_rate_calculations += length(i_cluster) + 1
            if is_IA_feasible(network, Partition(new_partition_matrix))
                # Fix BS j to this cluster
                F[:,j] = -Inf

                if length(i_cluster) == 1
                    # We have effectively put BS i in this cluster as well,
                    # without it really being aware. Do not try to put BS i
                    # in another cluster.
                    F[:,i] = -Inf
                end

                partition_matrix = new_partition_matrix
            else
                # Do not try to add BS j to the cluster of BS i again
                F[i_cluster,j] = -Inf
            end
        end
    end
    partition = Partition(partition_matrix)
    utilities, alphas, _ = longterm_utilities(channel, network, partition)
    a = restricted_growth_string(partition_matrix)
    objective = sum(utilities)
    Lumberjack.info("GreedyClustering finished.",
        @compat Dict(
            :sum_utility => objective,
            :a => a))

    # Store alphas as user priorities for precoding, if desired
    if aux_params["apply_overhead_prelog"]
        set_user_priorities!(network, alphas)
    end

    # Store cluster assignment together with existing cell assignment
    network.assignment = Assignment(temp_cell_assignment.cell_assignment, cluster_assignment_matrix(network, partition))

    # Return results
    results = AssignmentResults()
    results["utilities"] = utilities
    results["a"] = a
    results["alphas"] = alphas
    results["num_clusters"] = 1 + maximum(a)
    results["num_utility_calculations"] = num_utility_calculations
    results["num_longterm_rate_calculations"] = num_longterm_rate_calculations
    return results
end
