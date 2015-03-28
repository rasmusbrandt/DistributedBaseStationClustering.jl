# Performs branch and bound on a tree whose leaves describe restricted
# growth strings. The utilities are taken from utilities.jl, as well as
# the utopian bounds.
function BranchAndBoundClustering(channel, network)
    I = get_no_BSs(network)

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
    live = initialize_live(utopian_utilities); no_iters = 0; no_utility_calculations = 0
    while length(live) > 0
        no_iters += 1

        # Select next node to be processed. We use the best first strategy,
        # i.e. we pick the live node with the lowest (best) upper bound.
        sort!(live)
        parent = shift!(live)

        # Store incumbent evolution per iteration
        push!(incumbent_sum_utility_evolution, incumbent_sum_utility)

        for child in branch(parent)
            bound!(child, channel, network, utopian_utilities, I, assignment)
            no_utility_calculations += 1

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
    Lumberjack.info("BranchAndBoundClustering finished.",
        { :sum_utility => incumbent_sum_utility,
          :a => incumbent_a,
          :no_iters => no_iters,
          :no_utility_calculations => no_utility_calculations }
    )

    # Store cluster assignment together with existing cell assignment
    temp_cell_assignment = get_assignment(network)
    network.assignment = Assignment(temp_cell_assignment.cell_assignment, cluster_assignment_matrix(network, Partition(incumbent_a)))

    # Return results
    results = AssignmentResults()
    results["utilities"] = longterm_utilities(channel, network, Partition(incumbent_a))[1]
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
function bound!(node, channel, network, utopian_utilities, I, assignment)
    # The partial cluster is given by a
    partial_partition = Partition(node.a, skip_check=true) # By construction, a is a valid restricted growth string.

    # Rates for MSs already in clusters. These are utility bounds, since
    # the out-of-cluster interference of the unclustered users are not
    # taken into account.
    utility_bounds, _ = longterm_utilities(channel, network, partial_partition)

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
