##########################################################################
# Optimal base station clustering based on branch and bound.
#
# A branch and bound tree is developed which describes all possible
# restricted growth strings that describe all possible base station
# clusters.

function BranchAndBoundClustering(channel, network)
    # Perform cell selection
    LargeScaleFadingCellAssignment!(channel, network)
    assignment = get_assignment(network)

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
    apply_overhead_prelog = aux_params["apply_overhead_prelog"]
    IA_infeasible_negative_inf_utility = aux_params["IA_infeasible_negative_inf_utility"]

    # Consistency checks
    if I >= 12
        Lumberjack.warn("BranchAndBoundClustering will be slow since I = $I.")
    end
    if aux_params["clustering_type"] != :spectrum_sharing
        Lumberjack.error("BranchAndBoundClustering only works with spectrum sharing clustering.")
    end

    # Lumberjack.debug("BranchAndBoundClustering started.")

    # Utility lower bounds by using GreedyClustering_Multiple as initial incumbent.
    greedy_results = GreedyClustering_Multiple(channel, network)
    incumbent_utilities = greedy_results["utilities"]
    incumbent_a = greedy_results["a"]
    incumbent_sum_utility = sum(incumbent_utilities)
    # Lumberjack.debug("Incumbent (greedy) utilities calculated.", { :incumbent_utilities => incumbent_utilities, :incumbent_sum_utility => incumbent_sum_utility })

    # Perform eager branch and bound
    incumbent_sum_utility_evolution = Float64[]
    live = initialize_live(channel, network, Ps, sigma2s, I, Kc, M, N, d, assignment, apply_overhead_prelog, IA_infeasible_negative_inf_utility, E1_bound_in_rate_bound)
    no_iters = 0; no_utility_calculations = 0; no_longterm_rate_calculations = 0
    while length(live) > 0
        no_iters += 1

        # Select next node to be processed. We use the best first strategy,
        # i.e. we pick the live node with the highest (best) upper bound.
        parent = Base.Collections.heappop!(live, Base.Order.Reverse)

        # Store incumbent evolution per iteration
        push!(incumbent_sum_utility_evolution, incumbent_sum_utility)

        for child in branch(parent)
            bound!(child, channel, network, Ps, sigma2s, I, Kc, M, N, d, assignment, apply_overhead_prelog, IA_infeasible_negative_inf_utility, E1_bound_in_rate_bound)
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
function initialize_live(channel, network, Ps, sigma2s, I, Kc, M, N, d, assignment, apply_overhead_prelog, IA_infeasible_negative_inf_utility, E1_bound_in_rate_bound)
    root = BranchAndBoundNode([0], Inf)
    bound!(root, channel, network, Ps, sigma2s, I, Kc, M, N, d, assignment, apply_overhead_prelog, IA_infeasible_negative_inf_utility, E1_bound_in_rate_bound)
    # Lumberjack.debug("Root created.", { :node => root })
    return [ root ]
end

# Bound works by optimistically removing interference for unclustered BSs.
function bound!(node, channel, network, Ps, sigma2s, I, Kc, M, N, d, assignment, apply_overhead_prelog, IA_infeasible_negative_inf_utility, E1_bound_in_rate_bound)
    # The BSs are grouped based on their position in the graph. If they are put
    # in a cluster, they are 'clustered', otherwise they are 'unclustered'.
    all_BSs = IntSet(1:I)
    N_clustered = length(node.a); N_unclustered = I - N_clustered
    clustered_BSs = IntSet(1:N_clustered); nonclustered_BSs = IntSet(N_clustered+1:I)

    # For looping over clusters, we create a pseudo partition, where the
    # unclustered BSs belong to singleton blocks.
    m = 1 + maximum(node.a)
    pseudo_partition = Partition(cat(1, node.a, m:(m + (I-N_clustered-1))), skip_check=true)

    # The pre-log factor is upper bounded by the current pseudo_partition,
    # since each BS contributes the least towards CSI acquisition when it is
    # in a singleton cluster. This is the case for the 'unclustered'
    # BSs considered here. For the BSs that are 'clustered', we
    # know exactly their current CSI acquisition contribution.
    if apply_overhead_prelog
        alpha = spectrum_sharing_prelog_factor(network, pseudo_partition)
    else
        alpha = 1.
    end

    # The number of IA slots available will be important in the bound.
    # This number is given by the closed form in Liu2013. We also store
    # the BSs that are in the 'IA full' clusters, i.e. clusters that cannot
    # accept any more BSs. This is also used in the bound.
    N_available_IA_slots = Dict{Block,Int64}()
    BSs_in_full_clusters = IntSet()
    max_cluster_size = int((M + N - d)/(Kc*d))
    for block in pseudo_partition.blocks
        N_available_ = max_cluster_size - length(block.elements)
        N_available_IA_slots[block] = N_available_
        if N_available_ <= 0
            union!(BSs_in_full_clusters, block.elements)
        end
    end

    # Bound the utilities
    node_is_leaf = is_leaf(node, I)
    utility_bounds = zeros(Float64, I*Kc, d)
    for block in pseudo_partition.blocks
        N_available_IA_slots_ = N_available_IA_slots[block]

        # If this cluster is overloaded, the utilities suffer significantly.
        if N_available_IA_slots_ < 0
            # This cluster is overloaded since it violates the IA feasibility.
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
            # This cluster is not overloaded.
            for i in block.elements; for k in served_MS_ids(i, assignment)
                desired_power = channel.large_scale_fading_factor[k,i]*channel.large_scale_fading_factor[k,i]*(Ps[i]/(Kc*d)) # don't use ^2 for performance reasons

                # Leaves get true utility, other nodes get bound.
                if node_is_leaf
                    # All BSs outside my cluster contribute irreducible interference.
                    irreducible_interference_power = 0.
                    for j in setdiff(all_BSs, block.elements)
                        irreducible_interference_power += channel.large_scale_fading_factor[k,j]*channel.large_scale_fading_factor[k,j]*Ps[j]
                    end
                    rho = desired_power/(sigma2s[k] + irreducible_interference_power)
                else
                    # The SNR bound depends on if this BS is clustered or not.
                    if i <= N_clustered
                        # This BS is clustered.

                        # Clustered BSs outside my cluster contribute irreducible interference.
                        irreducible_interference_power = 0.
                        for j in setdiff(clustered_BSs, block.elements)
                            irreducible_interference_power += channel.large_scale_fading_factor[k,j]*channel.large_scale_fading_factor[k,j]*Ps[j]
                        end

                        # The bound is now due to picking the N_available_IA_slots_ strongest interferers (that are not clustered)
                        # and assuming that this interference is reducible. This is a bound since we are not ensuring the
                        # disjointness of the clusters here.
                        reducible_interference_levels = Float64[]
                        for j in nonclustered_BSs
                            push!(reducible_interference_levels, channel.large_scale_fading_factor[k,j]*channel.large_scale_fading_factor[k,j]*Ps[j])
                        end
                        sort!(reducible_interference_levels, rev=true)
                        rho = desired_power/(sigma2s[k] + irreducible_interference_power + sum(reducible_interference_levels[N_available_IA_slots_+1:end]))
                    else
                        # This BS is not clustered.

                        # I cannot joint any BS that belong to a full cluster, so
                        # those BSs contribute irreducible interference.
                        irreducible_interference_power = 0.
                        for j in BSs_in_full_clusters
                            irreducible_interference_power += channel.large_scale_fading_factor[k,j]*channel.large_scale_fading_factor[k,j]*Ps[j]
                        end

                        # We now pick the N_available_IA_slots_ strongest interferers (which do not belong to full clusters),
                        # and assume that this interference is reducible.
                        reducible_interference_levels = Float64[]
                        for j in setdiff(setdiff(all_BSs, BSs_in_full_clusters), block.elements)
                            push!(reducible_interference_levels, channel.large_scale_fading_factor[k,j]*channel.large_scale_fading_factor[k,j]*Ps[j])
                        end
                        sort!(reducible_interference_levels, rev=true)
                        rho = desired_power/(sigma2s[k] + irreducible_interference_power + sum(reducible_interference_levels[N_available_IA_slots_+1:end]))
                    end
                end

                # Finally, we can also bound the calculation of
                # exp(1/rho)*E1(1/rho) if that is desired.
                if E1_bound_in_rate_bound && !node_is_leaf
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
