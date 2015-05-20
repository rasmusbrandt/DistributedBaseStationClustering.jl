##########################################################################
# Optimal base station clustering based on branch and bound.
#
# A branch and bound tree is developed which describes all possible
# restricted growth strings that describe all possible base station
# clusters.

function BranchAndBoundClustering(channel, network)
    # We require symmetric networks, so we can write down the IA feasibility
    # condition in closed-form in the bound.
    check_Liu2013_applicability(network)
    check_Liu2013_symmetry(network)

    # Network parameters in symmetric network
    I = get_no_BSs(network); K = get_no_MSs(network)
    Kc = int(K/I)
    M = get_no_MS_antennas(network)[1]; N = get_no_BS_antennas(network)[1]
    d = get_no_streams(network)[1]
    Ps = get_transmit_powers(network); sigma2s = get_receiver_noise_powers(network)

    aux_params = get_aux_assignment_params(network)
    @defaultize_param! aux_params "BranchAndBoundClustering:E1_bound_in_rate_bound" false
    E1_bound_in_rate_bound = aux_params["BranchAndBoundClustering:E1_bound_in_rate_bound"]
    IA_infeasible_negative_inf_utility = aux_params["IA_infeasible_negative_inf_utility"]

    # Consistency checks
    if I >= 12
        Lumberjack.warn("BranchAndBoundClustering will be slow since I = $I.")
    end
    if aux_params["clustering_type"] != :spectrum_sharing
        Lumberjack.error("BranchAndBoundClustering only works with spectrum sharing clustering.")
    end

    # Lumberjack.debug("BranchAndBoundClustering started.")

    # Perform cell selection
    LargeScaleFadingCellAssignment!(channel, network)
    assignment = get_assignment(network)

    # Utility lower bounds by using GreedyClustering_Multiple as initial incumbent.
    greedy_results = GreedyClustering_Multiple(channel, network)
    incumbent_utilities = greedy_results["utilities"]
    incumbent_a = greedy_results["a"]
    incumbent_sum_utility = sum(incumbent_utilities)
    # Lumberjack.debug("Incumbent (greedy) utilities calculated.", { :incumbent_utilities => incumbent_utilities, :incumbent_sum_utility => incumbent_sum_utility })

    # Perform eager branch and bound
    incumbent_sum_utility_evolution = Float64[]
    live = initialize_live(channel, network, Ps, sigma2s, I, Kc, M, N, d, assignment, IA_infeasible_negative_inf_utility, E1_bound_in_rate_bound)
    no_iters = 0; no_utility_calculations = 0; no_longterm_rate_calculations = 0
    while length(live) > 0
        no_iters += 1

        # Select next node to be processed. We use the best first strategy,
        # i.e. we pick the live node with the highest (best) upper bound.
        parent = Base.Collections.heappop!(live, Base.Order.Reverse)

        # Store incumbent evolution per iteration
        push!(incumbent_sum_utility_evolution, incumbent_sum_utility)

        for child in branch(parent)
            bound!(child, channel, network, Ps, sigma2s, I, Kc, M, N, d, assignment, IA_infeasible_negative_inf_utility, E1_bound_in_rate_bound)
            no_utility_calculations += K
            no_longterm_rate_calculations += sum(child.a .== child.a[end]) # number of BSs affected by the child joining cluster a[end]

            # Is it worthwhile investigating this subtree/leaf more?
            if child.upper_bound > incumbent_sum_utility
                if is_leaf(child, I)
                    # For leaves, the upper bound is tight. Thus, we
                    # have found a new incumbent!
                    incumbent_sum_utility = child.upper_bound
                    incumbent_a = child.a

                    # Lumberjack.debug("Found new incumbent solution.",
                    #     { :node => child, :incumbent_sum_utility => incumbent_sum_utility }
                    # )
                else
                    # Lumberjack.debug("Keeping node since upper bound is above incumbent value.",
                    #     { :node => child, :incumbent_sum_utility => incumbent_sum_utility }
                    # )
                    Base.Collections.heappush!(live, child, Base.Order.Reverse)
                end
            else
                # Lumberjack.debug("Discarding node since upper bound is below incumbent value.",
                #     { :node => child, :incumbent_sum_utility => incumbent_sum_utility }
                # )
            end
        end
    end

    # Calculate final alphas
    final_partition = Partition(incumbent_a)
    utilities, alphas, _ = longterm_utilities(channel, network, final_partition)

    Lumberjack.info("BranchAndBoundClustering finished.",
        { :sum_utility => incumbent_sum_utility,
          :no_evaluated_partitions => no_utility_calculations/K,
          :a => incumbent_a }
    )

    # Store alphas as user priorities for precoding, if desired
    if aux_params["apply_overhead_prelog"]
        set_user_priorities!(network, alphas)
    end

    # Store cluster assignment together with existing cell assignment
    temp_cell_assignment = get_assignment(network)
    network.assignment = Assignment(temp_cell_assignment.cell_assignment, cluster_assignment_matrix(network, final_partition))

    # Return results
    results = AssignmentResults()
    results["utilities"] = utilities
    results["a"] = incumbent_a
    results["alphas"] = alphas
    results["no_clusters"] = 1 + maximum(incumbent_a)
    results["no_iters"] = no_iters
    results["no_utility_calculations"] = no_utility_calculations
    results["no_longterm_rate_calculations"] = no_longterm_rate_calculations
    return results
end

type BranchAndBoundNode
    a::Vector{Int} # restricted growth string describing the (partial) partition
    upper_bound::Float64
end

# For sorting of the list of live nodes
Base.isless(N1::BranchAndBoundNode, N2::BranchAndBoundNode) = (N1.upper_bound < N2.upper_bound)

# Helper functions
is_leaf(node, I) = (length(node.a) == I)

# Initialize the live structure by creating the root node.
function initialize_live(channel, network, Ps, sigma2s, I, Kc, M, N, d, assignment, IA_infeasible_negative_inf_utility, E1_bound_in_rate_bound)
    root = BranchAndBoundNode([0], Inf)
    bound!(root, channel, network, Ps, sigma2s, I, Kc, M, N, d, assignment, IA_infeasible_negative_inf_utility, E1_bound_in_rate_bound)
    # Lumberjack.debug("Root created.", { :node => root })
    return [ root ]
end

# Bound a node by testing feasibility and calculating the utilities for the
# clustered BSs and unclustered BSs.
function bound!(node, channel, network, Ps, sigma2s, I, Kc, M, N, d, assignment, IA_infeasible_negative_inf_utility, E1_bound_in_rate_bound)
    utility_bounds = zeros(Float64, I*Kc, d)

    # Create pseudo-cluster where the unclustered BSs are in singletons
    N_already_clustered = length(node.a); m = 1 + maximum(node.a)
    partition = Partition(cat(1, node.a, m:(m + (I-N_already_clustered-1))), skip_check=true) # valid restricted growth string by construction.

    # Calculate throughput bounds for all MSs
    for block in partition.blocks # N.B. all BSs are included!
        # The bound is based on the number of available 'IA slots' left
        # in this cluster. This is given by Liu's closed form.
        N_free_IA_slots = int((M + N - d)/(Kc*d) - length(block))

        if N_free_IA_slots < 0
            # IA infeasible partition
            if IA_infeasible_negative_inf_utility
                for i in block.elements; for k in served_MS_ids(i, assignment)
                    utility_bounds[k,:] = -Inf
                end; end
            else
                for i in block.elements; for k in served_MS_ids(i, assignment)
                    utility_bounds[k,:] = 0
                end; end
            end
        else
            # IA feasible partition
            intercluster_interferers = setdiff(IntSet(1:I), block.elements)
            N_intercluster_interferers = length(intercluster_interferers)
            for i in block.elements; for k in served_MS_ids(i, assignment)
                desired_power = channel.large_scale_fading_factor[k,i]*channel.large_scale_fading_factor[k,i]*(Ps[i]/(Kc*d)) # don't user ^2 for performance reasons

                # Bound the SNR if we are not at a leaf
                if is_leaf(node, I)
                    # Sum all interference
                    intercluster_interference_levels = Float64[]
                    for j in intercluster_interferers
                        Base.Collections.heappush!(intercluster_interference_levels, channel.large_scale_fading_factor[k,j]*channel.large_scale_fading_factor[k,j]*Ps[j], Base.Order.Reverse)
                    end
                    rho = desired_power/(sigma2s[k] + sum(intercluster_interference_levels))
                else
                    # Local SNR bound for this MS. This is a bound since we are not
                    # directly checking the disjoint condition required in the
                    # partition.
                    if N_free_IA_slots >= N_intercluster_interferers
                        # We can accommodate all interferers in this cluster.
                        rho = desired_power/sigma2s[k]
                    else
                        # We need to pick the N_free_IA_slots strongest interferers to include.
                        intercluster_interference_levels = Float64[]
                        for j in intercluster_interferers
                            Base.Collections.heappush!(intercluster_interference_levels, channel.large_scale_fading_factor[k,j]*channel.large_scale_fading_factor[k,j]*Ps[j], Base.Order.Reverse)
                        end
                        rho = desired_power/(sigma2s[k] + sum(intercluster_interference_levels[N_free_IA_slots+1:end]))
                    end
                end

                # alpha is upper bounded by the current partition, since each
                # BS contributes the least towards CSI acquisition when it is
                # in a singleton cluster. This is the case for the 'unclustered'
                # BSs considered here. For the BSs that are 'clustered', we
                # know exactly their current CSI acquisition contribution.
                alpha = spectrum_sharing_prelog_factor(network, partition)

                # Finally, we can also bound the calculation of
                # exp(1/rho)*E1(1/rho) if that is wanted.
                if E1_bound_in_rate_bound && !is_leaf(node, I)
                    utility_bounds[k,:] = alpha*longterm_rate(rho, :upper)
                else
                    utility_bounds[k,:] = alpha*longterm_rate(rho, :none)
                end
            end; end
        end
    end

    # The final sum utility bound is obtained by summing the individual bounds.
    node.upper_bound = sum(utility_bounds)

    # Lumberjack.debug("Bounding.", { :node => node })
end

# Branch a node by creating a number of descendants, corresponding to putting
# the BS at this depth in different clusters. Branched nodes inherit the
# parent bound, until their bounds are updated.
function branch(parent)
    m = 1 + maximum(parent.a)

    no_children = m + 1
    children = Array(BranchAndBoundNode, no_children)

    for p = 1:no_children
        child_a = push!(copy(parent.a), p - 1)
        child = BranchAndBoundNode(child_a, parent.upper_bound)
        children[p] = child

        # Lumberjack.debug("Branching.", { :node => child })
    end

    return children
end
