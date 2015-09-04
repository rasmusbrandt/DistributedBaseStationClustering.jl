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
    I = get_num_BSs(network); K = get_num_MSs(network)
    Kc = int(K/I)
    M = get_num_BS_antennas(network)[1]; N = get_num_MS_antennas(network)[1]
    d = get_num_streams(network)[1]
    Ps = get_transmit_powers(network); sigma2s = get_receiver_noise_powers(network)

    aux_params = get_aux_assignment_params(network)
    num_coherence_symbols = get_aux_network_param(network, "num_coherence_symbols")
    beta_network_sdma = get_aux_network_param(network, "beta_network_sdma")
    @defaultize_param! aux_params "BranchAndBoundClustering:max_abs_optimality_gap" 0.
    max_abs_optimality_gap = aux_params["BranchAndBoundClustering:max_abs_optimality_gap"]
    @defaultize_param! aux_params "BranchAndBoundClustering:max_rel_optimality_gap" 0.
    max_rel_optimality_gap = aux_params["BranchAndBoundClustering:max_rel_optimality_gap"]
    @defaultize_param! aux_params "BranchAndBoundClustering:E1_bound_in_rate_bound" false
    E1_bound_in_rate_bound = aux_params["BranchAndBoundClustering:E1_bound_in_rate_bound"]
    @defaultize_param! aux_params "BranchAndBoundClustering:store_fathomed_subtree_sizes" false
    store_fathomed_subtree_sizes = aux_params["BranchAndBoundClustering:store_fathomed_subtree_sizes"]

    # Lumberjack.debug("BranchAndBoundClustering started.")

    # Throughput lower bounds by trying different heuristics.
    incumbent_throughputs = zeros(Float64, K, d)
    incumbent_a = zeros(Int, I)
    incumbent_sum_throughput = -Inf
    for heuristic in [NoClustering, GrandCoalitionClustering, GreedyClustering_Single, GreedyClustering_Multiple]
        heuristic_results = heuristic(channel, network)
        heuristic_sum_throughput = sum(heuristic_results["throughputs"])
        if heuristic_sum_throughput > incumbent_sum_throughput
            incumbent_throughputs = heuristic_results["throughputs"]
            incumbent_a = heuristic_results["a"]
            incumbent_sum_throughput = heuristic_sum_throughput
        end
        # Lumberjack.debug("Potential incumbent throughputs calculated.", { :heuristic => heuristic, :incumbent_sum_throughput => incumbent_sum_throughput })
    end

    # Optimal solution from previous branch and bound round
    if has_aux_network_param(network, "BranchAndBoundClustering:cache:optimal_a")
        previous_a = get_aux_network_param(network, "BranchAndBoundClustering:cache:optimal_a")
        previous_partition = Partition(previous_a)
        previous_throughputs, = longterm_throughputs(channel, network, previous_partition) # with current power allocation
        previous_sum_throughput = sum(previous_throughputs)
        if previous_sum_throughput > incumbent_sum_throughput
            incumbent_throughputs = previous_throughputs
            incumbent_a = previous_a
            incumbent_sum_throughput = previous_sum_throughput
        end
    end

    # Perform eager branch and bound
    lower_bound_evolution = Float64[]; upper_bound_evolution = Float64[]; fathoming_evolution = Int[]
    live = initialize_live(channel, network, Ps, sigma2s, I, Kc, M, N, d, beta_network_sdma, num_coherence_symbols, assignment, E1_bound_in_rate_bound)
    num_iters = 0; num_sum_throughput_calculations = 0
    abs_conv_crit = 0.; premature_ending = false
    while length(live) > 0
        num_iters += 1

        # Select next node to be processed. We use the best first strategy,
        # i.e. we pick the live node with the highest (best) upper bound.
        parent = Base.Collections.heappop!(live, Base.Order.Reverse)

        # Store bound evolution per iteration
        push!(lower_bound_evolution, incumbent_sum_throughput)
        push!(upper_bound_evolution, parent.upper_bound)

        # Check convergence (parent has the highest upper bound)
        abs_conv_crit = parent.upper_bound - incumbent_sum_throughput
        rel_conv_crit = abs_conv_crit/incumbent_sum_throughput
        if abs_conv_crit <= max_abs_optimality_gap || rel_conv_crit <= max_rel_optimality_gap
            # Lumberjack.debug("Converged.", { :abs_conv_crit => abs_conv_crit, :max_abs_optimality_gap => max_abs_optimality_gap })
            premature_ending = true
            break
        end

        fathomed_subtree_size = 0
        for child in branch(parent)
            bound!(child, channel, network, Ps, sigma2s, I, Kc, M, N, d, beta_network_sdma, num_coherence_symbols, assignment, E1_bound_in_rate_bound)
            num_sum_throughput_calculations += 1

            # Is it worthwhile investigating this subtree/leaf more?
            if child.upper_bound > incumbent_sum_throughput
                if is_leaf(child, I)
                    # For leaves, the upper bound is tight. Thus, we
                    # have found a new incumbent!
                    incumbent_sum_throughput = child.upper_bound
                    incumbent_a = child.a

                    # Lumberjack.debug("Found new incumbent solution.",
                    #     { :node => child, :incumbent_sum_throughput => incumbent_sum_throughput }
                    # )
                else
                    # Lumberjack.debug("Keeping node since upper bound is above incumbent value.",
                    #     { :node => child, :incumbent_sum_throughput => incumbent_sum_throughput }
                    # )
                    Base.Collections.heappush!(live, child, Base.Order.Reverse)
                end
            else
                # Lumberjack.debug("Discarding node since upper bound is below incumbent value.",
                #     { :node => child, :incumbent_sum_throughput => incumbent_sum_throughput }
                # )

                if store_fathomed_subtree_sizes
                    fathomed_subtree_size += subtree_size(child, I)
                end
            end
        end

        if store_fathomed_subtree_sizes
            push!(fathoming_evolution, fathomed_subtree_size)
        end
    end

    # Did we find the global optimum?
    if (abs_conv_crit < 0.) || !premature_ending
        abs_conv_crit = 0.
    end

    # Calculate final prelogs
    final_partition = Partition(incumbent_a)
    throughputs, throughputs_split, _, prelogs = longterm_throughputs(channel, network, final_partition)

    # Verify that our local throughput calculation makes sens
    sum_throughput = sum(throughputs)
    if abs(sum_throughput - incumbent_sum_throughput) > 1e-10
        Lumberjack.error("Something is wrong in the bound.")
    end

    Lumberjack.info("BranchAndBoundClustering finished.",
        { :sum_throughput => sum_throughput,
          :num_sum_throughput_calculations => num_sum_throughput_calculations,
          :abs_conv_crit => abs_conv_crit,
          :max_abs_optimality_gap => max_abs_optimality_gap,
          :a => incumbent_a }
    )

    # Store a for next branch and bound run
    set_aux_network_param!(network, incumbent_a, "BranchAndBoundClustering:cache:optimal_a")

    # Store prelogs for precoding
    set_aux_network_param!(network, prelogs[1], "prelogs_cluster_sdma")
    set_aux_network_param!(network, prelogs[2], "prelogs_network_sdma")

    # Store cluster assignment together with existing cell assignment
    temp_cell_assignment = get_assignment(network)
    network.assignment = Assignment(temp_cell_assignment.cell_assignment, cluster_assignment_matrix(network, final_partition))

    # Return results
    results = AssignmentResults()
    results["throughputs"] = throughputs
    results["throughputs_cluster_sdma"] = throughputs_split[1]
    results["throughputs_network_sdma"] = throughputs_split[2]
    results["a"] = incumbent_a
    results["num_clusters"] = 1 + maximum(incumbent_a)
    results["avg_cluster_size"] = avg_cluster_size(incumbent_a)
    results["num_sum_throughput_calculations"] = num_sum_throughput_calculations
    results["num_iters"] = num_iters
    results["lower_bound_evolution"] = reshape(lower_bound_evolution, (1, 1, length(lower_bound_evolution)))
    results["upper_bound_evolution"] = reshape(upper_bound_evolution, (1, 1, length(upper_bound_evolution)))
    results["fathoming_evolution"] = reshape(fathoming_evolution, (1, 1, length(fathoming_evolution)))
    return results
end

type BranchAndBoundNode
    a::Vector{Int} # restricted growth string describing the (partial) partition
    upper_bound::Float64
end

# For sorting of the list of live nodes
Base.isless(N1::BranchAndBoundNode, N2::BranchAndBoundNode) = (N1.upper_bound < N2.upper_bound)

# Helper functions
m(node) = 1 + maximum(node.a)
depth(node) = size(node.a, 1)
is_leaf(node, I) = (depth(node) == I)
num_children(node) = 1 + m(node)

subtree_size(node, I) = subtree_size(depth(node), m(node), I)
function subtree_size(depth, m, I)
    if depth == I
        return 1
    else
        return 1 + m*subtree_size(depth+1, m, I) + subtree_size(depth+1, m+1 ,I)
    end
end

# Initialize the live structure by creating the root node.
function initialize_live(channel, network, Ps, sigma2s, I, Kc, M, N, d, beta_network_sdma, num_coherence_symbols, assignment, E1_bound_in_rate_bound)
    root = BranchAndBoundNode([0], Inf)
    bound!(root, channel, network, Ps, sigma2s, I, Kc, M, N, d, beta_network_sdma, num_coherence_symbols, assignment, E1_bound_in_rate_bound)
    # Lumberjack.debug("Root created.", { :node => root })
    return [ root ]
end

# Bound works by optimistically removing interference for unclustered BSs.
function bound!(node, channel, network, Ps, sigma2s, I, Kc, M, N, d,
    beta_network_sdma, num_coherence_symbols, assignment, E1_bound_in_rate_bound)

    beta_cluster_sdma = 1 - beta_network_sdma

    # The BSs are grouped based on their position in the graph. If they are put
    # in a cluster, they are 'clustered', otherwise they are 'unclustered'.
    all_BSs = IntSet(1:I)
    N_clustered = length(node.a); N_unclustered = I - N_clustered
    clustered_BSs = IntSet(1:N_clustered); nonclustered_BSs = IntSet(N_clustered+1:I)

    # For looping over clusters, we create a pseudo partition, where the
    # unclustered BSs belong to singleton blocks.
    pseudo_partition_a = Array(Int64, I)
    for i = 1:N_clustered
        pseudo_partition_a[i] = node.a[i]
    end
    m = 1 + maximum(node.a)
    for i = (N_clustered+1):I
        pseudo_partition_a[i] = m + (i - N_clustered - 1)
    end
    pseudo_partition = Partition(pseudo_partition_a, skip_check=true)

    # The number of IA slots available will be important in the bounds.
    # This number is given by the closed form in Liu2013. We also store
    # the BSs that are in the 'IA full' clusters, i.e. clusters that cannot
    # accept any more BSs. This is also used in the bound.
    max_cluster_size = ifloor((M + N - d)/(Kc*d))
    N_available_IA_slots = Dict{Block,Int64}()
    BSs_in_full_clusters = IntSet()
    for block in pseudo_partition.blocks
        N_available_ = max_cluster_size - length(block.elements)
        N_available_IA_slots[block] = N_available_
        if N_available_ <= 0
            union!(BSs_in_full_clusters, block.elements)
        end
    end
    BSs_in_nonfull_clusters = setdiff(all_BSs, BSs_in_full_clusters)
    N_BSs_in_nonfull_clusters = length(BSs_in_nonfull_clusters)

    # Find cluster SDMA prelog optimal cluster size (to be used in prelog bounds)
    num_symbols_cluster_sdma = beta_cluster_sdma*num_coherence_symbols
    cs_opt = min((num_symbols_cluster_sdma - I*(M + Kc*(N + d)))/(2*I*Kc*M), max_cluster_size)
    cs1 = iceil(cs_opt); cs2 = ifloor(cs_opt) # cs_opt might be fractional so check surrounding integers. The overhead function is unimodal, it's OK.
    a1 = symmetric_prelog_cluster_sdma(cs1, beta_cluster_sdma, num_symbols_cluster_sdma, I, M, Kc, N, d)
    a2 = symmetric_prelog_cluster_sdma(cs2, beta_cluster_sdma, num_symbols_cluster_sdma, I, M, Kc, N, d)
    optimal_cluster_size_cluster_sdma = (a1 > a2) ? cs1 : cs2

    # Prelog, rate, and throughput bounds
    node_is_leaf = is_leaf(node, I)
    prelog_bounds_cluster_sdma = zeros(Float64, I*Kc)
    prelog_bounds_network_sdma = zeros(Float64, I*Kc)
    throughput_bounds = zeros(Float64, I*Kc, d)
    for block in pseudo_partition.blocks
        N_available_IA_slots_ = N_available_IA_slots[block]

        # Clusters that are IA overloaded get zero throughput, both for cluster SDMA and for network SDMA.
        if N_available_IA_slots_ >= 0
            cluster_size = length(block.elements)

            # For the prelog bound, we want the number of BSs in this cluster
            # to be close to the optimal number.
            if cluster_size < optimal_cluster_size_cluster_sdma
                # There is still space, so add more BSs to this cluster. We never want more than optimal_cluster_size_cluster_sdma
                # BSs in our clusters, because that is when the prelog starts going down again.
                clustered_BS_cluster_size_bound = min(cluster_size + min(N_available_IA_slots_, N_unclustered), optimal_cluster_size_cluster_sdma)
                nonclustered_BS_cluster_size_bound = min(cluster_size + min(N_available_IA_slots_, N_BSs_in_nonfull_clusters), optimal_cluster_size_cluster_sdma)
            else
                # We (potentially) already have too many BSs in this cluster,
                # thus our prelog can only go down by adding more. Bound the
                # prelog by our current value.
                clustered_BS_cluster_size_bound = cluster_size
                nonclustered_BS_cluster_size_bound = cluster_size
            end
            clustered_BS_prelog_bound_cluster_sdma = symmetric_prelog_cluster_sdma(clustered_BS_cluster_size_bound, beta_cluster_sdma, num_symbols_cluster_sdma, I, M, Kc, N, d)
            nonclustered_BS_prelog_bound_cluster_sdma = symmetric_prelog_cluster_sdma(nonclustered_BS_cluster_size_bound, beta_cluster_sdma, num_symbols_cluster_sdma, I, M, Kc, N, d)

            # True prelog when all BSs are clustered.
            leaf_prelog_cluster_sdma = symmetric_prelog_cluster_sdma(cluster_size, beta_cluster_sdma, num_symbols_cluster_sdma, I, M, Kc, N, d)

            # Stuff needed for rate bounds
            outside_all_BSs = setdiff(all_BSs, block.elements)
            outside_clustered_BSs = setdiff(clustered_BSs, block.elements)
            outside_BSs_in_nonfull_clusters = setdiff(BSs_in_nonfull_clusters, block.elements)

            # Get appropriate bounds for all MSs in this cluster.
            for i in block.elements; for k in served_MS_ids(i, assignment)
                # Desired channel.
                desired_power = channel.large_scale_fading_factor[k,i]*channel.large_scale_fading_factor[k,i]*(Ps[i]/(Kc*d)) # don't use ^2 for performance reasons
                rho_cluster_sdma = desired_power/sigma2s[k]

                # Leaves get true values, other nodes get bound.
                if node_is_leaf
                    # True prelog.
                    prelog_bounds_cluster_sdma[k] = leaf_prelog_cluster_sdma

                    # All BSs outside my cluster contribute irreducible interference.
                    irreducible_interference_power = 0.
                    for j in outside_all_BSs
                        irreducible_interference_power += channel.large_scale_fading_factor[k,j]*channel.large_scale_fading_factor[k,j]*Ps[j]
                    end
                    rho_network_sdma = desired_power/(sigma2s[k] + irreducible_interference_power)
                else
                    # The bounds depends on if this BS is clustered or not.
                    if i <= N_clustered
                        # This BS is clustered.

                        # Bound cluster SDMA prelog by adding suitably many BSs to this cluster
                        prelog_bounds_cluster_sdma[k] = clustered_BS_prelog_bound_cluster_sdma

                        # Clustered BSs outside my cluster contribute irreducible interference.
                        irreducible_interference_power = 0.
                        for j in outside_clustered_BSs
                            irreducible_interference_power += channel.large_scale_fading_factor[k,j]*channel.large_scale_fading_factor[k,j]*Ps[j]
                        end

                        # The bound is now due to picking the N_available_IA_slots_ strongest interferers (that are not clustered)
                        # and assuming that this interference is reducible. This is a bound since we are not ensuring the
                        # disjointness of the clusters here.
                        reducible_interference_levels = Float64[]
                        for j in nonclustered_BSs
                            push!(reducible_interference_levels, channel.large_scale_fading_factor[k,j]*channel.large_scale_fading_factor[k,j]*Ps[j])
                        end
                        sort!(reducible_interference_levels, rev=true) # Could speed this up by using a heap.
                        rho_network_sdma = desired_power/(sigma2s[k] + irreducible_interference_power + sum(reducible_interference_levels[N_available_IA_slots_+1:end]))
                    else
                        # This BS is not clustered.

                        # Bound cluster SDMA prelog by adding suitably many BSs to this cluster
                        prelog_bounds_cluster_sdma[k] = nonclustered_BS_prelog_bound_cluster_sdma

                        # I cannot joint any BS that belong to a full cluster, so
                        # those BSs contribute irreducible interference.
                        irreducible_interference_power = 0.
                        for j in BSs_in_full_clusters
                            irreducible_interference_power += channel.large_scale_fading_factor[k,j]*channel.large_scale_fading_factor[k,j]*Ps[j]
                        end

                        # We now pick the N_available_IA_slots_ strongest interferers (which do not belong to full clusters),
                        # and assume that this interference is reducible.
                        reducible_interference_levels = Float64[]
                        for j in outside_BSs_in_nonfull_clusters
                            push!(reducible_interference_levels, channel.large_scale_fading_factor[k,j]*channel.large_scale_fading_factor[k,j]*Ps[j])
                        end
                        sort!(reducible_interference_levels, rev=true) # Could speed this up by using a heap.
                        rho_network_sdma = desired_power/(sigma2s[k] + irreducible_interference_power + sum(reducible_interference_levels[N_available_IA_slots_+1:end]))
                    end
                end

                # Network SDMA is either zero (CSI acquisition infeasible)
                # or fixed to the beta constant.
                if prelog_bounds_cluster_sdma[k] == 0.
                    prelog_bounds_network_sdma[k] = 0.
                else
                    prelog_bounds_network_sdma[k] = beta_network_sdma
                end

                # Finally, we can also bound the calculation of
                # exp(1/rho)*E1(1/rho) if that is desired. (Note that the
                # prelogs are bounded above.)
                if E1_bound_in_rate_bound && !node_is_leaf
                    throughput_bounds[k,:] =
                        prelog_bounds_cluster_sdma[k]*exp_times_E1(rho_cluster_sdma, bound=:upper) +
                        prelog_bounds_network_sdma[k]*exp_times_E1(rho_network_sdma, bound=:upper)
                else
                    throughput_bounds[k,:] =
                        prelog_bounds_cluster_sdma[k]*exp_times_E1(rho_cluster_sdma) +
                        prelog_bounds_network_sdma[k]*exp_times_E1(rho_network_sdma)
                end
            end; end
        end
    end

    # The final sum throughput bound is obtained by summing the individual bounds.
    node.upper_bound = sum(throughput_bounds)

    # Lumberjack.debug("Bounding.", { :node => node })
end

# Branch a node by creating a number of descendants, corresponding to putting
# the BS at this depth in different clusters. Branched nodes inherit the
# parent bound, until their bounds are updated.
function branch(parent)
    num_children_ = num_children(parent)
    children = Array(BranchAndBoundNode, num_children_)

    for p = 1:num_children_
        child_a = push!(copy(parent.a), p - 1)
        child = BranchAndBoundNode(child_a, parent.upper_bound)
        children[p] = child

        # Lumberjack.debug("Branching.", { :node => child })
    end

    return children
end

# Special case of the cluster SDMA prelog used a lot in the bound.
symmetric_prelog_cluster_sdma(cluster_size, beta_cluster_sdma, num_symbols_cluster_sdma, I, M, Kc, N, d) =
    beta_cluster_sdma*max(0., cluster_size/I - (cluster_size*(M + Kc*(N + d)) + cluster_size^2*Kc*M)/num_symbols_cluster_sdma)
