##########################################################################
# Optimal base station clustering based on branch and bound.
#
# A branch and bound tree is developed which describes all possible
# restricted growth strings that describe all possible base station
# clusters. The utilities are taken from utilities.jl. The bounds are
# very simple: The pre-log factors are bounded by the current pre-log
# factor for already clustered BSs, and bounded by 1 for non-clustered BSs.
# The spectral efficiency bounds are assuming that non-clustered BSs do
# not contribute interference to already clustered BSs, and non-clustered
# BSs have their utopian (interference-free) spectral efficiencies as bounds.

function BranchAndBoundClustering(channel, network)
    I = get_num_BSs(network); K = get_num_MSs(network)

    aux_params = get_aux_assignment_params(network)
    @defaultize_param! aux_params "BranchAndBoundClustering:bracket_E1" false
    bracket_E1 = aux_params["BranchAndBoundClustering:bracket_E1"]

    # Warn if this will be slow...
    if I >= 12
        Lumberjack.warn("BranchAndBoundClustering will be slow since I = $I.")
    end

    # Lumberjack.debug("BranchAndBoundClustering started.")

    # Perform cell selection
    LargeScaleFadingCellAssignment!(channel, network)
    assignment = get_assignment(network)

    # Utility upper bounds by assuming feasibility of grand coalition.
    grand_coalition_a = zeros(Int, I) # restricted growth string with all zeros
    grand_coalition = Partition(grand_coalition_a)
    _, _, utopian_utilities = longterm_utilities(channel, network, grand_coalition)
    # utopian_sum_value = sum(utopian_utilities)
    # Lumberjack.debug("Utopian (fully cooperative) utilities calculated.", { :utopian_utilities => utopian_utilities, :utopian_sum_value => utopian_sum_value })

    # Utility lower bounds by using GreedyClustering as initial incumbent.
    greedy_results = GreedyClustering(channel, network)
    incumbent_utilities = greedy_results["utilities"]
    incumbent_a = greedy_results["a"]
    incumbent_sum_utility = sum(incumbent_utilities)
    # Lumberjack.debug("Incumbent (greedy) utilities calculated.", { :incumbent_utilities => incumbent_utilities, :incumbent_sum_utility => incumbent_sum_utility })

    # Perform eager branch and bound
    incumbent_sum_utility_evolution = Float64[]
    live = initialize_live(utopian_utilities);
    num_iters = 0
    num_utility_calculations = 0; num_longterm_rate_calculations = 0
    while length(live) > 0
        num_iters += 1

        # Select next node to be processed. We use the best first strategy,
        # i.e. we pick the live node with the highest (best) upper bound.
        sort!(live, rev=true)
        parent = shift!(live)

        # Store incumbent evolution per iteration
        push!(incumbent_sum_utility_evolution, incumbent_sum_utility)

        for child in branch(parent)
            bound!(child, channel, network, utopian_utilities, I, assignment, bracket_E1)
            num_utility_calculations += K
            num_longterm_rate_calculations += sum(child.a .== child.a[end]) # number of BSs affected by the child joining cluster a[end]

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
                    push!(live, child)
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
        @compat Dict(
            :sum_utility => incumbent_sum_utility,
            :a => incumbent_a))

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
    results["num_clusters"] = 1 + maximum(incumbent_a)
    results["num_iters"] = num_iters
    results["num_utility_calculations"] = num_utility_calculations
    results["num_longterm_rate_calculations"] = num_longterm_rate_calculations
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
function initialize_live(utopian_utilities)
    root = BranchAndBoundNode([0], sum(utopian_utilities))
    # Lumberjack.debug("Root created.", { :node => root })
    return [ root ]
end

# Bound a node by testing feasibility and calculating the utilities for the
# clustered BSs and unclustered BSs.
function bound!(node, channel, network, utopian_utilities, I, assignment, bracket_E1)
    # The partial cluster is given by a
    partial_partition = Partition(node.a, skip_check=true) # By construction, a is a valid restricted growth string.

    # Rates for MSs already in clusters. These are utility bounds, since
    # the out-of-cluster interference of the unclustered users are not
    # taken into account.
    if is_leaf(node, I) || !bracket_E1
        utility_bounds, _ = longterm_utilities(channel, network, partial_partition)
    else
        utility_bounds, _ = longterm_utilities(channel, network, partial_partition, bound=:upper)
    end

    # Bound the unclustered users utilities by their utopian utilities.
    for j in setdiff(IntSet(1:I), IntSet(1:length(node.a))); for l in served_MS_ids(j, assignment)
        utility_bounds[l,:] = utopian_utilities[l,:]
    end; end

    # Sum utility bound. This might be -Inf, if any of the utilities are -Inf.
    # This can happen when a particular block is not IA feasible, and
    # the aux_assignment_param IA_infeasible_utility_inf is set to true.
    node.upper_bound = sum(utility_bounds)

    # Lumberjack.debug("Bounding.", { :node => node })
end

# Branch a node by creating a number of descendants, corresponding to putting
# the BS at this depth in different clusters. Branched nodes inherit the
# parent bound, until their bounds are updated.
function branch(parent)
    m = 1 + maximum(parent.a)

    num_children = m + 1
    children = Array(BranchAndBoundNode, num_children)

    for p = 1:num_children
        child_a = push!(copy(parent.a), p - 1)
        child = BranchAndBoundNode(child_a, parent.upper_bound)
        children[p] = child

        # Lumberjack.debug("Branching.", { :node => child })
    end

    return children
end
