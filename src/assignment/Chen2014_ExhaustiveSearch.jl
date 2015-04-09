##########################################################################
# Optimal base station clustering based on exhaustive search over the
# utilities proposed in the paper
#
# Chen, Cheng, "Clustering for Interference Alignment in Multiuser
# Interference Network," IEEE Trans. Vehicular Technology, vol. 63, no. 6,
# pp. 2613-2624, July 2014, doi: 10.1109/TVT.2013.2292897

function Chen2014_ExhaustiveSearch(channel, network)
    I = get_no_BSs(network); K = get_no_MSs(network)
    aux_params = get_aux_assignment_params(network)

    # Consistency check
    if I != K
        Lumberjack.error("Chen2014_LinearObjClustering can only handle I = K scenarios.")
    end
    if aux_params["IA_infeasible_negative_inf_utility"] == false
        Lumberjack.info("Chen2014_ExhaustiveSearch only finds solutions where all clusters are IA feasible. IA_infeasible_negative_inf_utility is set to false, which means that the other methods might find solutions where some blocks are turned off due to IA infeasibility.")
    end
    if I > 12
        Lumberjack.warn("Chen2014_ExhaustiveSearch will be slow since I = $I.")
    end

    # Perform cell selection
    LargeScaleFadingCellAssignment!(channel, network)

    # Get W matrix
    W = Chen2014_W_matrix(channel, network)

    # Exhaustive search over partitions
    no_iters = 0
    best_objective = 0.
    best_partition = Partition()
    for partition in PartitionIterator(I)
        no_iters += 1

        # Check that IA is feasible for this cluster structure. Note that this
        # means that Chen2014_ExhaustiveSearch cannot handle situations where
        # IA infeasible blocks are turned off, e.g. when the aux_assignment_param
        # IA_infeasible_negative_inf_utility is set to false.
        if is_IA_feasible(network, partition)
            # Calculate objective
            objective = 0.
            for block in partition.blocks
                objective = 0.
                for i in block.elements; for j in block.elements
                    objective += W[i,j]
                end; end
            end

            if objective > best_objective
                best_objective = objective
                best_partition = partition
            end
        end
    end
    utilities, alphas, _ = longterm_utilities(channel, network, best_partition)
    a = restricted_growth_string(best_partition)
    Lumberjack.info("Chen2014_ExhaustiveSearch finished.",
        { :sum_utility => sum(utilities),
          :a => a,
          :alphas => alphas,
          :no_iters => no_iters,
          :Chen2014_objective => best_objective }
    )

    # Store alphas as user priorities for precoding, if desired
    if aux_params["apply_overhead_prelog"]
        set_user_priorities!(network, alphas)
    end

    # Store cluster assignment together with existing cell assignment
    temp_cell_assignment = get_assignment(network)
    network.assignment = Assignment(network.assignment.cell_assignment, cluster_assignment_matrix(network, best_partition))

    # Return results
    results = AssignmentResults()
    results["utilities"] = utilities
    results["a"] = a
    results["alphas"] = alphas
    results["no_iters"] = no_iters
    results["no_clusters"] = 1 + maximum(a)
    results["Chen2014_objective"] = best_objective
    return results
end

function Chen2014_W_matrix(channel, network)
    I = get_no_BSs(network)
    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    assignment = get_assignment(network)

    W = zeros(Float64, I, I)
    for i = 1:I; for j = 1:I
        if i == j
            continue
        end

        k = collect(served_MS_ids(i, assignment))[1]
        desired_power = channel.large_scale_fading_factor[k,i]^2*Ps[i]*channel.large_scale_fading_factor[k,j]^2*Ps[j]
        int_noise_power = sigma2s[k] + channel.large_scale_fading_factor[k,j]^2*Ps[j] + channel.large_scale_fading_factor[k,i]^2*Ps[i]
        W[i,j] = log2(1 + desired_power/int_noise_power)

        l = collect(served_MS_ids(j, assignment))[1]
        desired_power = channel.large_scale_fading_factor[l,j]^2*Ps[j]*channel.large_scale_fading_factor[l,i]^2*Ps[i]
        int_noise_power = sigma2s[k] + channel.large_scale_fading_factor[l,i]^2*Ps[i] + channel.large_scale_fading_factor[l,j]^2*Ps[j]
        W[i,j] += log2(1 + desired_power/int_noise_power)
    end; end

    return W
end
