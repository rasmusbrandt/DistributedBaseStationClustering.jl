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
    Kc = convert(Int, K/I)
    M = get_num_BS_antennas(network)[1]; N = get_num_MS_antennas(network)[1]
    d = get_num_streams(network)[1]
    Ps = get_transmit_powers(network); sigma2s = get_receiver_noise_powers(network)

    aux_params = get_aux_assignment_params(network)
    num_coherence_symbols = get_aux_network_param(network, "num_coherence_symbols")
    beta_network_sdma = get_aux_network_param(network, "beta_network_sdma")
    @defaultize_param! aux_params "BranchAndBoundClustering:branching_rule" :dfs
    branching_rule = aux_params["BranchAndBoundClustering:branching_rule"]
    @defaultize_param! aux_params "BranchAndBoundClustering:max_abs_optimality_gap" 0.
    max_abs_optimality_gap = aux_params["BranchAndBoundClustering:max_abs_optimality_gap"]
    @defaultize_param! aux_params "BranchAndBoundClustering:max_rel_optimality_gap" 0.
    max_rel_optimality_gap = aux_params["BranchAndBoundClustering:max_rel_optimality_gap"]
    @defaultize_param! aux_params "BranchAndBoundClustering:E1_bound_in_rate_bound" false
    E1_bound_in_rate_bound = aux_params["BranchAndBoundClustering:E1_bound_in_rate_bound"]
    @defaultize_param! aux_params "BranchAndBoundClustering:store_evolution" false
    store_evolution = aux_params["BranchAndBoundClustering:store_evolution"]

    # Lumberjack.debug("BranchAndBoundClustering started.")

    # Best upper bound (mainly needed for DFS)
    best_upper_bound = Inf

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

    # Precompute desired powers and cluster SDMA spectral efficiencies
    desired_powers = Array(Float64, I*Kc)
    rates_cluster_sdma = Array(Float64, I*Kc)
    for i = 1:I; for k in served_MS_ids(i, assignment)
        desired_powers[k] = channel.large_scale_fading_factor[k,i]*channel.large_scale_fading_factor[k,i]*(Ps[i]/(Kc*d)) # don't use ^2 for performance reasons
        rho_cluster_sdma = desired_powers[k]/sigma2s[k]
        if E1_bound_in_rate_bound
            rates_cluster_sdma[k] = exp_times_E1(rho_cluster_sdma, bound=:upper)
        else
            rates_cluster_sdma[k] = exp_times_E1(rho_cluster_sdma)
        end
    end; end

    # Precompute interfering powers
    interfering_powers = Array(Float64, I, I*Kc)
    for i = 1:I; for k in served_MS_ids(i, assignment)
        for j = 1:I
            interfering_powers[j,k] = channel.large_scale_fading_factor[k,j]*channel.large_scale_fading_factor[k,j]*Ps[j]
        end
    end; end

    # Perform eager branch and bound
    lower_bound_evolution = Float64[]; upper_bound_evolution = Float64[]; fathoming_evolution = Int[]
    live = initialize_live(channel, network, Ps, sigma2s, I, Kc, M, N, d, beta_network_sdma, num_coherence_symbols, assignment, E1_bound_in_rate_bound, desired_powers, interfering_powers, rates_cluster_sdma)
    num_iters = 0; num_bounded_nodes = 0
    abs_conv_crit = 0.; premature_ending = false
    while length(live) > 0
        num_iters += 1

        # Select next node to be processed.
        if branching_rule == :bfs
            # Best first, i.e. the highest (best) upper bound
            parent = Base.Collections.heappop!(live, Base.Order.Reverse)
            best_upper_bound = parent.upper_bound
        elseif branching_rule == :dfs
            # Depth first.
            best_upper_bound = maximum([ node.upper_bound for node in live ])
            parent = pop!(live)
        end

        if store_evolution
            # Store bound evolution per iteration
            push!(lower_bound_evolution, incumbent_sum_throughput)
            push!(upper_bound_evolution, best_upper_bound)
        end

        # Check convergence
        abs_conv_crit = best_upper_bound - incumbent_sum_throughput
        rel_conv_crit = abs_conv_crit/incumbent_sum_throughput
        if abs_conv_crit <= max_abs_optimality_gap || rel_conv_crit <= max_rel_optimality_gap
            # Lumberjack.debug("Converged.", { :abs_conv_crit => abs_conv_crit, :max_abs_optimality_gap => max_abs_optimality_gap })
            premature_ending = true
            break
        end

        fathomed_subtree_size = 0
        for child in branch(parent)
            bound!(child, channel, network, Ps, sigma2s, I, Kc, M, N, d, beta_network_sdma, num_coherence_symbols, assignment, E1_bound_in_rate_bound, desired_powers, interfering_powers, rates_cluster_sdma)
            num_bounded_nodes += 1

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
                    if branching_rule == :bfs
                        Base.Collections.heappush!(live, child, Base.Order.Reverse)
                    else
                        push!(live, child)
                    end
                end
            else
                # Lumberjack.debug("Discarding node since upper bound is below incumbent value.",
                #     { :node => child, :incumbent_sum_throughput => incumbent_sum_throughput }
                # )

                if store_evolution
                    fathomed_subtree_size += subtree_size(child, I) - 1 # minus one since we already explored child
                end
            end
        end

        if store_evolution
            push!(fathoming_evolution, fathomed_subtree_size)
        end
    end

    # Add remaining subtrees that were implicitly fathomed
    if store_evolution
        push!(fathoming_evolution, sum([ subtree_size(node, I) for node in live ]))
    end

    # Did we find the global optimum?
    if (abs_conv_crit < 0.) || !premature_ending
        abs_conv_crit = 0.
    end

    # Calculate final prelogs
    final_partition = Partition(incumbent_a)
    throughputs, throughputs_split, _, prelogs = longterm_throughputs(channel, network, final_partition)

    # Verify that our local throughput calculation makes sense
    sum_throughput = sum(throughputs)
    if abs(sum_throughput - incumbent_sum_throughput) > 1e-10
        Lumberjack.error("Something is wrong in the bound.")
    end

    Lumberjack.info("BranchAndBoundClustering finished.",
        @Compat.Dict(
            :sum_throughput => sum_throughput,
            :num_bounded_nodes => num_bounded_nodes,
            :abs_conv_crit => abs_conv_crit,
            :max_abs_optimality_gap => max_abs_optimality_gap,
            :a => incumbent_a)
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
    results["num_clusters"] = reshape(1 + maximum(incumbent_a), 1, 1)
    results["avg_cluster_size"] = reshape(avg_cluster_size(incumbent_a), 1, 1)
    results["num_iters"] = reshape(num_iters, 1, 1)
    results["num_bounded_nodes"] = reshape(num_bounded_nodes, 1, 1)
    results["num_sum_throughput_calculations"] = reshape(num_bounded_nodes, 1, 1)
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
function initialize_live(channel, network, Ps, sigma2s, I, Kc, M, N, d,
    beta_network_sdma, num_coherence_symbols, assignment,
    E1_bound_in_rate_bound, desired_powers, interfering_powers, rates_cluster_sdma)

    root = BranchAndBoundNode([0], Inf)
    bound!(root, channel, network, Ps, sigma2s, I, Kc, M, N, d, beta_network_sdma, num_coherence_symbols, assignment, E1_bound_in_rate_bound, desired_powers, interfering_powers, rates_cluster_sdma)
    # Lumberjack.debug("Root created.", { :node => root })
    return [ root ]
end

# Bound works by optimistically removing interference for unclustered BSs.
function bound!(node, channel, network, Ps, sigma2s, I, Kc, M, N, d,
    beta_network_sdma, num_coherence_symbols, assignment,
    E1_bound_in_rate_bound, desired_powers, interfering_powers, rates_cluster_sdma)

    # Preallocate some memory
    scratch = Array(Float64, I)
    aggregated = fill(0., 2)

    beta_cluster_sdma = 1 - beta_network_sdma

    # The BSs are grouped based on their position in the graph. If they are put
    # in a cluster, they are 'clustered', otherwise they are 'unclustered'.
    all_BSs = IntSet(1:I)
    N_clustered = length(node.a); N_unclustered = I - N_clustered
    clustered_BSs = IntSet(1:N_clustered); unclustered_BSs = IntSet(N_clustered+1:I)
    reducible_interference_levels1 = Array(Float64, N_unclustered)

    # For looping over clusters, we create a pseudo partition, where the
    # unclustered BSs belong to singleton blocks.
    pseudo_partition_a = Array(Int, I)
    @simd for i = 1:N_clustered
        @inbounds pseudo_partition_a[i] = node.a[i]
    end
    m = 1 + maximum(node.a)
    @simd for i = (N_clustered+1):I
        @inbounds pseudo_partition_a[i] = m + (i - N_clustered - 1)
    end
    pseudo_partition = Partition(pseudo_partition_a, skip_check=true)

    # The number of IA slots available will be important in the bounds.
    # This number is given by the closed form in Liu2013. We also store
    # the BSs that are in the 'IA full' clusters, i.e. clusters that cannot
    # accept any more BSs. This is also used in the bound.
    max_cluster_size = floor(Int, (M + N - d)/(Kc*d))
    BSs_in_full_clusters = IntSet()
    for block in pseudo_partition.blocks
        N_available_ = max_cluster_size - length(block.elements)
        if N_available_ <= 0
            union!(BSs_in_full_clusters, block.elements)
        end
    end
    BSs_in_nonfull_clusters = setdiff(all_BSs, BSs_in_full_clusters)
    N_BSs_in_nonfull_clusters = length(BSs_in_nonfull_clusters)

    # Find cluster SDMA prelog optimal cluster size (to be used in prelog bounds)
    num_symbols_cluster_sdma = beta_cluster_sdma*num_coherence_symbols
    cs_opt = min((num_symbols_cluster_sdma - I*(M + Kc*(N + d)))/(2*I*Kc*M), max_cluster_size)
    cs1 = ceil(Int, cs_opt); cs2 = floor(Int, cs_opt) # cs_opt might be fractional so check surrounding integers. The overhead function is unimodal, it's OK.
    a1 = symmetric_prelog_cluster_sdma(cs1, beta_cluster_sdma, num_symbols_cluster_sdma, I, M, Kc, N, d)
    a2 = symmetric_prelog_cluster_sdma(cs2, beta_cluster_sdma, num_symbols_cluster_sdma, I, M, Kc, N, d)
    optimal_cluster_size_cluster_sdma = (a1 > a2) ? cs1 : cs2

    # Prelog, rate, and throughput bounds
    node_is_leaf = is_leaf(node, I)
    prelog_bounds_cluster_sdma = zeros(Float64, I*Kc)
    prelog_bounds_network_sdma = zeros(Float64, I*Kc)
    throughput_bounds = zeros(Float64, I*Kc, d)
    for block in pseudo_partition.blocks
        cluster_size = length(block.elements)
        N_available_IA_slots_ = max_cluster_size - cluster_size

        # Clusters that are IA overloaded get zero throughput, both for cluster SDMA and for network SDMA.
        if N_available_IA_slots_ >= 0
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
            N_outside_BSs_in_nonfull_clusters = length(outside_BSs_in_nonfull_clusters)

            # Get appropriate bounds for all MSs in this cluster.
            for i in block.elements; for k in served_MS_ids(i, assignment)
                @inbounds aggregated[1] = 0.; @inbounds aggregated[2] = 0.

                # Leaves get true values, other nodes get bound.
                if node_is_leaf
                    # True prelog.
                    prelog_bounds_cluster_sdma[k] = leaf_prelog_cluster_sdma

                    # All BSs outside my cluster contribute irreducible interference.
                    sum_irreducible_interference_power!(aggregated, k, outside_all_BSs, interfering_powers)
                else
                    # The bounds depends on if this BS is clustered or not.
                    if i <= N_clustered
                        # This BS is clustered.

                        # Bound cluster SDMA prelog by adding suitably many BSs to this cluster
                        prelog_bounds_cluster_sdma[k] = clustered_BS_prelog_bound_cluster_sdma

                        # Clustered BSs outside my cluster contribute irreducible interference.
                        sum_irreducible_interference_power!(aggregated, k, outside_clustered_BSs, interfering_powers)

                        # The bound is now due to picking the N_available_IA_slots_ strongest interferers (that are not clustered)
                        # and assuming that this interference is reducible. This is a bound since we are not ensuring the
                        # disjointness of the clusters here.
                        sum_reducible_interference_power!(aggregated, k, unclustered_BSs, N_unclustered - N_available_IA_slots_,
                                                          interfering_powers, sub(scratch, 1:N_unclustered))
                    else
                        # This BS is not clustered.

                        # Bound cluster SDMA prelog by adding suitably many BSs to this cluster
                        prelog_bounds_cluster_sdma[k] = nonclustered_BS_prelog_bound_cluster_sdma

                        # I cannot joint any BS that belong to a full cluster, so
                        # those BSs contribute irreducible interference.
                        sum_irreducible_interference_power!(aggregated, k, BSs_in_full_clusters, interfering_powers)

                        # We now pick the N_available_IA_slots_ strongest interferers (which do not belong to full clusters),
                        # and assume that this interference is reducible.
                        sum_reducible_interference_power!(aggregated, k, outside_BSs_in_nonfull_clusters, N_outside_BSs_in_nonfull_clusters - N_available_IA_slots_,
                                                          interfering_powers, sub(scratch, 1:N_outside_BSs_in_nonfull_clusters))
                    end
                end
                @inbounds rho_network_sdma = desired_powers[k]/(sigma2s[k] + sum(aggregated))

                # Network SDMA is either zero (CSI acquisition infeasible)
                # or fixed to the beta constant.
                @inbounds begin
                    if prelog_bounds_cluster_sdma[k] == 0.
                        prelog_bounds_network_sdma[k] = 0.
                    else
                        prelog_bounds_network_sdma[k] = beta_network_sdma
                    end
                end

                # Finally, we can also bound the calculation of
                # exp(1/rho)*E1(1/rho) if that is desired. (Note that the
                # prelogs are bounded above.)
                if E1_bound_in_rate_bound && !node_is_leaf
                    @inbounds throughput_bound =
                        prelog_bounds_cluster_sdma[k]*rates_cluster_sdma[k] + # rates_cluster_sdma is already bounded in the calling function
                        prelog_bounds_network_sdma[k]*exp_times_E1(rho_network_sdma, bound=:upper)
                    @simd for n = 1:d
                        @inbounds throughput_bounds[k,n] = throughput_bound
                    end
                else
                    throughput_bound =
                        prelog_bounds_cluster_sdma[k]*rates_cluster_sdma[k] +
                        prelog_bounds_network_sdma[k]*exp_times_E1(rho_network_sdma)
                    @simd for n = 1:d
                        @inbounds throughput_bounds[k,n] = throughput_bound
                    end
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

function sum_irreducible_interference_power!(aggregated, k, from_BSs, interfering_powers)
    for j in from_BSs
        @inbounds aggregated[1] += interfering_powers[j,k]
    end
end

function sum_reducible_interference_power!(aggregated, k, from_BSs, N_available, interfering_powers, scratch)
    idx = 1
    for j in from_BSs
        @inbounds scratch[idx] = interfering_powers[j,k]
        idx += 1
    end
    sort!(scratch)
    @simd for idx in 1:N_available
        @inbounds aggregated[2] += scratch[idx]
    end
end

# Special case of the cluster SDMA prelog used a lot in the bound.
symmetric_prelog_cluster_sdma(cluster_size, beta_cluster_sdma, num_symbols_cluster_sdma, I, M, Kc, N, d) =
    beta_cluster_sdma*max(0., cluster_size/I - (cluster_size*(M + Kc*(N + d)) + cluster_size^2*Kc*M)/num_symbols_cluster_sdma)
