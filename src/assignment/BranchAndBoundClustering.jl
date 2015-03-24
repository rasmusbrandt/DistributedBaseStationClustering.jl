type BranchAndBoundNode
    a::Vector{Int} # restricted growth string describing the (partial) partition
    upper_bound::Float64
end

# For sorting of the list of live nodes
Base.isless(N1::BranchAndBoundNode, N2::BranchAndBoundNode) = (N1.upper_bound < N2.upper_bound)

# Helper functions
is_leaf(node, I) = (length(node.a) == I)

# Initialize the live structure by creating the root node.
function initialize_live(utopian_rates)
    root = BranchAndBoundNode([0], sum(utopian_rates))
    Lumberjack.debug("Root created.", { :node => root })
    return [ root ]
end

# Bound a node by testing feasibility and calculating the rates for the
# clustered BSs and unclustered BSs.
function bound!(node, channel, network, utopian_rates)
    I = get_no_BSs(network)
    assignment = get_assignment(network)

    # The partial cluster is given by a
    partial_partition = Partition(node.a)

    if is_IA_feasible(network, partial_partition)
        # This is a feasible allocation. Bound rates.

        # Rates for MSs already in clusters. These are rate bounds, since
        # the out-of-cluster interference of the unclustered users are not
        # taken into account.
        rate_bounds = longterm_throughputs(channel, network, partial_partition)

        # Bound the unclustered users rates by their utopian rates.
        for j in setdiff(1:I, 1:length(node.a)); for l in served_MS_ids(j, assignment)
            rate_bounds[l] = utopian_rates[l]
        end; end

        # Sum rate bound
        node.upper_bound = sum(rate_bounds)
    else
        # Infeasible allocation. 
        node.upper_bound = -Inf
    end

    Lumberjack.debug("Bounding.", { :node => node })
end

# Branch a node by creating a number of descendants, corresponding to putting
# the BS at this depth in different clusters. Branched nodes inherit the
# parent bound, until their bounds are updated.
function branch(parent)
    m = 1 + maximum(parent.a)

    no_children = m + 1
    children = Array(BranchAndBoundNode, no_children)

    for p = 1:no_children
        child_a = vcat(parent.a, p - 1)
        child = BranchAndBoundNode(child_a, parent.upper_bound)
        children[p] = child

        Lumberjack.debug("Branching.", { :node => child })
    end

    return children
end

function BranchAndBoundClustering(channel, network)
    I = get_no_BSs(network)

    Lumberjack.debug("BranchAndBoundClustering started.")

    # Perform cell selection
    LargeScaleFadingCellAssignment!(channel, network)

    # Rate upper bounds by assuming feasibility of grand coalition.
    grand_coalition_a = zeros(Int, I) # restricted growth string with all zeros
    grand_coalition = Partition(grand_coalition_a)
    utopian_rates = longterm_throughputs(channel, network, grand_coalition)
    utopian_value = sum(utopian_rates)
    Lumberjack.debug("Utopian (fully cooperative) rates calculated.", { :utopian_rates => utopian_rates, :utopian_value => utopian_value })

    # Incumbent: non-cooperative case
    incumbent_a, incumbent_rates = GreedyClustering(channel, network)
    incumbent_value = sum(incumbent_rates)
    Lumberjack.debug("Incumbent (greedy) rates calculated.", { :incumbent_rates => incumbent_rates, :incumbent_value => incumbent_value })

    # Perform eager branch and bound
    incumbent_evolution = Float64[]
    live = initialize_live(utopian_rates); no_iters = 0; no_feasibility_checks = 0
    while length(live) > 0
        no_iters += 1

        # Select next node to be processed. (Nodes with the least upper bounds
        # should be investigated first.)
        sort!(live)
        parent = shift!(live)

        # Store incumbent evolution per investigated node
        push!(incumbent_evolution, incumbent_value)

        for child in branch(parent)
            bound!(child, channel, network, utopian_rates)
            no_feasibility_checks += 1

            # Is it worthwhile investigating this subtree/leaf more?
            if child.upper_bound > incumbent_value
                if is_leaf(child, I)
                    # For leaves, the upper bound is tight. Thus, we
                    # have found a new incumbent!
                    incumbent_value = child.upper_bound
                    incumbent_a = child.a

                    Lumberjack.debug("Found new incumbent solution.",
                        { :node => child, :incumbent_value => incumbent_value }
                    )
                else
                    Lumberjack.debug("Keeping node since upper bound is above incumbent value.",
                        { :node => child, :incumbent_value => incumbent_value }
                    )
                    push!(live, child)
                end
            else
                Lumberjack.debug("Discarding node since upper bound is below incumbent value.",
                    { :node => child, :incumbent_value => incumbent_value }
                )
            end
        end
    end
    Lumberjack.info("BranchAndBoundClustering finished.",
        { :sum_rate => incumbent_value,
          :a => incumbent_a,
          :no_iters => no_iters,
          :no_feasibility_checks => no_feasibility_checks }
    )

    # Store cluster assignment together with existing cell assignment
    temp_cell_assignment = get_assignment(network)
    network.assignment = Assignment(temp_cell_assignment.cell_assignment, cluster_assignment_matrix(network, Partition(incumbent_a)))
end
