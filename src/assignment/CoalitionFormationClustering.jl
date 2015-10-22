##########################################################################
# Base station clustering based on coalitional formation. Two solution
# concepts are used: individual-based stability and group-based stability.
#
# In both algorithms, we have *BS_utilities*, which are the MS utilities
# for the served MSs summed up. (This *could* be generalized to other functions.)
#
# In the individual-based stability algorithm (CoalitionFormationClustering_Individual),
# each BS is allowed to deviate to another coalition, given that the BSs
# in that coalition accepts the deviating BS. Letting all BSs consecutively
# deviate eventually leads to a coalition structure which is individually stable.
# The parameter that this algorithm takes are:
#
#   Maximum number of deviation searches that each BS is allowed to performed
#       CoalitionFormationClustering_Individual:search_budget (Int)
#
#   The way the possible deviations are ordered
#       CoalitionFormationClustering_Individual:search_order
#   This can be either :greedy or :fair, where :greedy means that the BS
#   with the largest utility is allowed to deviate first, and :fair means
#   that the BS with the smallest utility is allowed to deviate first.
#
# The group-based stability algorithm (CoalitionFormationClustering_Group) works
# directly on coalitions, rather than on BSs as the individual-based stability
# algorithm did. In this case, coalitions are allowed to merge, leading up to
# a coalition structure which is group stable.
# The parameter that this algorithm takes are:
#
#   Maximum number of coalitions merging
#       CoalitionFormationClustering_Group:max_num_merging_coalitions (Int)
#
#   The way the possible mergers are ordered
#       CoalitionFormationClustering_Group:search_order
#   This can be either :greedy or :lexicographic, where :greedy means that
#   the merger which leads to the largest sum utility (over the entire network)
#   will be tried first. :lexicographic means that the mergers are tried
#   in lexicographic order.

##########################################################################
# Individual-based stability algorithm
type CoalitionFormationClustering_IndividualState
    partition::Partition
    BS_utilities::Vector{Float64}
    history::Vector{Set{IntSet}}
    num_searches::Vector{Int}
    num_utility_calculations::Int
    num_longterm_rate_calculations::Int
    movie_state::MovieState
end

function CoalitionFormationClustering_Individual(channel, network)
    I = get_num_BSs(network); K = get_num_MSs(network)

    aux_params = get_aux_assignment_params(network)
    @defaultize_param! aux_params "CoalitionFormationClustering_Individual:search_budget" 10
    @defaultize_param! aux_params "CoalitionFormationClustering_Individual:search_order" :lexicographic
    @defaultize_param! aux_params "CoalitionFormationClustering_Individual:stability_type" :individual
    @defaultize_param! aux_params "CoalitionFormationClustering_Individual:use_history" true
    @defaultize_param! aux_params "CoalitionFormationClustering_Individual:starting_point" :grand
    @defaultize_param! aux_params "CoalitionFormationClustering_Individual:prepare_movie" false
    search_budget = aux_params["CoalitionFormationClustering_Individual:search_budget"]
    search_order = aux_params["CoalitionFormationClustering_Individual:search_order"]
    stability_type = aux_params["CoalitionFormationClustering_Individual:stability_type"]
    use_history = aux_params["CoalitionFormationClustering_Individual:use_history"]
    starting_point = aux_params["CoalitionFormationClustering_Individual:starting_point"]
    prepare_movie = aux_params["CoalitionFormationClustering_Individual:prepare_movie"]

    # Perform cell selection
    LargeScaleFadingCellAssignment!(channel, network)
    temp_cell_assignment = get_assignment(network)

    # Initial coalition structure
    if starting_point == :grand
        initial_partition = Partition(zeros(Int, I))
    elseif starting_point == :singletons
        initial_partition = Partition(collect(0:(I-1)))
    end
    initial_BS_utilities = longterm_BS_utilities(channel, network, initial_partition, temp_cell_assignment, I)
    initial_num_searches = zeros(Int, I)
    state = CoalitionFormationClustering_IndividualState(
                initial_partition,
                initial_BS_utilities,
                [ Set{IntSet}() for i = 1:I ],
                initial_num_searches,
                K,
                K,
                MovieState())
    if prepare_movie
        push_event!(state.movie_state, :accept, "Initial coalition structure.", 0, [], [], initial_partition, initial_BS_utilities)
    end

    # Let each BS deviate, and stop when no BS deviates (individual-based stability)
    deviation_performed = trues(I) # temporary, to enter the loop
    while any(deviation_performed)
        # Give all BSs a chance to deviate. If search_order == :greedy, we
        # let the BSs deviate in the order of their current utilities, i.e. the
        # BS doing the best is going first. If search_order == :fair instead,
        # the BS which is doing the worst will go first. This is more similar
        # to GreedyClustering, where the strongest interfering links are clustered
        # first.
        if search_order == :greedy
            ordered_BS_list = sortperm(state.BS_utilities, rev=true)
        elseif search_order == :fair
            ordered_BS_list = sortperm(state.BS_utilities, rev=false)
        elseif search_order == :lexicographic
            ordered_BS_list = 1:I
        end

        deviation_performed = falses(I)
        for i in ordered_BS_list
            deviation_performed[i] = deviate!(state, i, I, K, search_budget,
                stability_type, use_history, prepare_movie, channel, network, temp_cell_assignment)
        end
    end
    push_event!(state.movie_state, :accepted_by_all, "Individual stability achieved.", 0, [], [], state.partition, state.BS_utilities)
    utilities, alphas, _ = longterm_utilities(channel, network, state.partition)
    a = restricted_growth_string(state.partition)
    Lumberjack.info("CoalitionFormationClustering_Individual finished.",
        @compat Dict(
            :sum_utility => sum(utilities),
            :a => a))

    # Store alphas as user priorities for precoding, if desired
    if aux_params["apply_overhead_prelog"]
        set_user_priorities!(network, alphas)
    end

    # Store cluster assignment together with existing cell assignment
    network.assignment = Assignment(temp_cell_assignment.cell_assignment, cluster_assignment_matrix(network, state.partition))

    # Return results
    results = AssignmentResults()
    results["utilities"] = utilities
    results["a"] = a
    results["alphas"] = alphas
    results["num_clusters"] = 1 + maximum(a)
    results["num_searches"] = state.num_searches
    results["num_utility_calculations"] = state.num_utility_calculations
    results["num_longterm_rate_calculations"] = state.num_longterm_rate_calculations
    if prepare_movie
        results["movie_state"] = state.movie_state
    end
    return results
end

# Lets BS i deviate in the individual stability model.
# Returns true if it did deviate, otherwise false.
function deviate!(state::CoalitionFormationClustering_IndividualState, i, I, K,
    search_budget, stability_type, use_history, prepare_movie,
    channel, network, cell_assignment)

    # First check that we have not exceeded our search budget
    if state.num_searches[i] >= search_budget
        return false
    end

    # Divide blocks such that old_block is the block that BS i used to
    # belong to, and other_blocks is an array of all other blocks.
    old_block = Block() # to store variable from inside loop
    other_blocks = Block[]
    for block in state.partition.blocks
        if i in block.elements
            old_block = Block(setdiff(block.elements, IntSet(i))) # BS i does not belong to the old_block anymore
        else
            push!(other_blocks, block)
        end
    end

    # Create all possible deviations for BS i
    BS_not_singleton_coalition_before = length(old_block) > 0 ? true : false # was this BS not in a singleton coalition before?
    num_new_partitions = length(other_blocks) + convert(Int, BS_not_singleton_coalition_before)
    new_partitions = Array(Partition, num_new_partitions)
    deviated_BS_utilities = zeros(Float64, I, num_new_partitions)
    for n = 1:length(other_blocks)
        # Loop over the deviations where BS i joins an existing coalition
        new_partition = Partition()
        other_blocks_cp = deepcopy(other_blocks) # need to deepcopy, so the created coalitions will not all be the same...
        union!(new_partition.blocks, other_blocks_cp)

        # Add the old block unless it used to be a singleton
        if BS_not_singleton_coalition_before
            push!(new_partition.blocks, old_block)
        end

        # Add BS i to coalition n
        push!(other_blocks_cp[n].elements, i)
        new_partitions[n] = new_partition
        deviated_BS_utilities[:,n] = longterm_BS_utilities(channel, network, new_partition, cell_assignment, I)

        # Complexity metrics
        state.num_utility_calculations += K
        state.num_longterm_rate_calculations += length(old_block) + length(other_blocks_cp[n])
    end
    if BS_not_singleton_coalition_before
        # BS i was in a non-singleton coalition before deviation. Add the the
        # possibility that it belongs to a non-singleton coalition after deviation.
        new_partition = Partition()
        other_blocks_cp = deepcopy(other_blocks)
        union!(new_partition.blocks, other_blocks_cp)

        # Add the old block
        push!(new_partition.blocks, old_block)

        # Add BS i to new singleton coalition
        push!(new_partition.blocks, Block(IntSet(i)))
        new_partitions[end] = new_partition
        deviated_BS_utilities[:,end] = longterm_BS_utilities(channel, network, new_partition, cell_assignment, I)

        # Complexity metrics
        state.num_utility_calculations += K
        state.num_longterm_rate_calculations += length(old_block) + 1
    end

    # Preliminary Nash stability check. No need to try to deviate unless
    # BS i improves in its utility.
    if !any(deviated_BS_utilities[i,:] .> state.BS_utilities[i])
        push_event!(state.movie_state, :accept, "Cell $i does not want to leave its current coalition.", i, [], [], state.partition, state.BS_utilities)
        return false
    end

    # Check deviations, trying to join the coalitions in the order that
    # benefits BS i the most.
    for sort_idx in sortperm(squeeze(deviated_BS_utilities[i,:], 1), rev=true)
        # Stop searching if we otherwise would exceed our search budget.
        if state.num_searches[i] >= search_budget
            return false
        end

        # Let's try to deviate
        state.num_searches[i] += 1

        # Find block that BS i belongs to in this partition
        my_block = Block()
        for block in new_partitions[sort_idx].blocks
            if i in block.elements
                my_block = block
                break
            end
        end

        # Check if the existing members of this coalition allow BS i to join
        # (this check includes BS i, unnecessarily)
        BSs_in_new_block = collect(my_block.elements); BSs_in_new_block_except_me = setdiff(BSs_in_new_block, IntSet(i))
        BSs_in_old_block = collect(old_block.elements)
        if prepare_movie
            # Nash criterion
            if deviated_BS_utilities[i,sort_idx] > state.BS_utilities[i]
                # History set
                if !use_history || !in(IntSet(BSs_in_new_block), state.history[i])
                    if length(BSs_in_new_block_except_me) == 0
                        push_event!(state.movie_state, :ask, "Cell $i wants to leave for its singleton coalition.", i, [], [], state.partition, state.BS_utilities)
                    else
                        push_event!(state.movie_state, :ask, "Cell $i wants to join coalition $BSs_in_new_block_except_me.", i, BSs_in_new_block_except_me, [], state.partition, state.BS_utilities)
                    end
                end
            end
        end
        if individual_stability(deviated_BS_utilities[:,sort_idx], state.BS_utilities, i, BSs_in_new_block, BSs_in_old_block, state.history[i], stability_type, use_history)
            # Let BS i join this coalition
            state.partition = new_partitions[sort_idx]
            state.BS_utilities = deviated_BS_utilities[:,sort_idx]

            # Add coalition to history
            push!(state.history[i], IntSet(BSs_in_new_block))

            if prepare_movie
                push_event!(state.movie_state, :accepted_by_all, "Cell $i was accepted by all.", i, BSs_in_new_block_except_me, trues(BSs_in_new_block_except_me), state.partition, state.BS_utilities)
            end

            return true
        elseif prepare_movie
            # Nash criterion
            if deviated_BS_utilities[i,sort_idx] > state.BS_utilities[i]
                # History set
                if !use_history || !in(IntSet(BSs_in_new_block), state.history[i])
                    answers = deviated_BS_utilities[BSs_in_new_block_except_me,sort_idx] .>= state.BS_utilities[BSs_in_new_block_except_me]
                    if all(answers .== false)
                        push_event!(state.movie_state, :rejected_by_all, "Cell $i was rejected by all.", i, BSs_in_new_block_except_me, answers, state.partition, state.BS_utilities)
                    else
                        push_event!(state.movie_state, :mixed_response, "Cell $i was rejected by $(BSs_in_new_block_except_me[!answers]) but accepted by $(BSs_in_new_block_except_me[answers]).", i, BSs_in_new_block_except_me, answers, state.partition, state.BS_utilities)
                    end
                end
            end
        end
    end
    return false
end

# Check stability of a particular deviating BS for the individual
# coalition formation algorithm.
function individual_stability(new_BS_utilities, old_BS_utilities,
    deviating_BS_idx, new_coalition_idxs, old_coalition_idxs, history,
    stability_type, use_history)

    # Check that this coalition does not exist in the history
    if use_history && IntSet(new_coalition_idxs) in history
        return false
    end

    # Check that BS i improves
    nash = (new_BS_utilities[deviating_BS_idx] > old_BS_utilities[deviating_BS_idx])
    if stability_type == :nash
        return nash
    end

    # Check if the BSs in the new coalition improve
    individual = (nash && all(new_BS_utilities[new_coalition_idxs] .>= old_BS_utilities[new_coalition_idxs]))
    if stability_type == :individual
        return individual
    end

    # Check if the BSs in the old coalition improve
    contractual = (individual && all(new_BS_utilities[old_coalition_idxs] .>= old_BS_utilities[old_coalition_idxs]))
    if stability_type == :contractual
        return contractual
    end
end

##########################################################################
# Group-based stability algorithm
type CoalitionFormationClustering_GroupState
    partition::Partition
    BS_utilities::Vector{Float64}
    r::Int
    num_utility_calculations::Int
    num_longterm_rate_calculations::Int
end

function CoalitionFormationClustering_Group(channel, network)
    I = get_num_BSs(network); K = get_num_MSs(network)

    aux_params = get_aux_assignment_params(network)
    @defaultize_param! aux_params "CoalitionFormationClustering_Group:max_num_merging_coalitions" 3
    @defaultize_param! aux_params "CoalitionFormationClustering_Group:search_order" :greedy
    max_num_merging_coalitions = aux_params["CoalitionFormationClustering_Group:max_num_merging_coalitions"]
    search_order = aux_params["CoalitionFormationClustering_Group:search_order"]
    in(search_order, [:greedy, :lexicographic]) || Lumberjack.error("Incorrect CoalitionFormationClustering_Group:search_order.")

    # Perform cell selection
    LargeScaleFadingCellAssignment!(channel, network)
    temp_cell_assignment = get_assignment(network)

    # Initial coalition structure is the non-cooperative state
    initial_partition = Partition(collect(0:(I-1)))
    initial_BS_utilities = longterm_BS_utilities(channel, network, initial_partition, temp_cell_assignment, I)
    state = CoalitionFormationClustering_GroupState(initial_partition, initial_BS_utilities, min(I, aux_params["CoalitionFormationClustering_Group:max_num_merging_coalitions"]), K, K)

    # Let coalitions merge, until no coalitions want to merge (group-based stability)
    while state.r >= 2 && length(state.partition) > 2
        # Keep merging until no coalitions want to merge
        merge_performed = true
        while merge_performed
            merge_performed = merge!(state, I, K, max_num_merging_coalitions,
                search_order, channel, network, temp_cell_assignment)
        end

        # No more mergers happened with the current r, so decrease it.
        state.r -= 1
    end
    utilities, alphas, _ = longterm_utilities(channel, network, state.partition)
    a = restricted_growth_string(state.partition)
    Lumberjack.info("CoalitionFormationClustering_Group finished.",
        @compat Dict(
            :sum_utility => sum(utilities),
            :a => a))

    # Store alphas as user priorities for precoding, if desired
    if aux_params["apply_overhead_prelog"]
        set_user_priorities!(network, alphas)
    end

    # Store cluster assignment together with existing cell assignment
    network.assignment = Assignment(temp_cell_assignment.cell_assignment, cluster_assignment_matrix(network, state.partition))

    # Return results
    results = AssignmentResults()
    results["utilities"] = utilities
    results["a"] = a
    results["alphas"] = alphas
    results["num_clusters"] = 1 + maximum(a)
    results["num_utility_calculations"] = state.num_utility_calculations
    results["num_longterm_rate_calculations"] = state.num_longterm_rate_calculations
    return results
end

# Enumerates all possible mergers, given the current r and coalition structure,
# and then tries to perform mergers in a specific order given by the params.
function merge!(state::CoalitionFormationClustering_GroupState, I, K,
    max_num_merging_coalitions, search_order, channel, network, cell_assignment)

    # Put current blocks in an array for easy indexing
    all_blocks = collect(state.partition.blocks)
    all_blocks_card = length(all_blocks)

    # Create all possible r-mergers
    num_new_partitions = binomial(all_blocks_card, state.r) # all_blocks_card choose r
    new_partitions = Array(Partition, num_new_partitions)
    merged_BS_utilities = zeros(Float64, I, num_new_partitions)
    merged_BSs = Array(Vector{Int}, num_new_partitions)

    # Loop lexicographically over all r-combinations of the current coalitions
    n = 0
    for merging_blocks_idxs in kCombinationIterator(all_blocks_card, state.r)
        n += 1

        # Separate blocks based on who is merging and not
        merged_block = Block()
        for merging_blocks_idx in merging_blocks_idxs
            union!(merged_block.elements, all_blocks[merging_blocks_idx].elements)
        end
        other_blocks = all_blocks[collect(setdiff(IntSet(1:all_blocks_card), IntSet(merging_blocks_idxs)))]

        # Partition after merging
        new_partition = Partition()
        push!(new_partition.blocks, merged_block)
        union!(new_partition.blocks, other_blocks)
        new_partitions[n] = new_partition
        merged_BSs[n] = collect(merged_block.elements)
        merged_BS_utilities[:,n] = longterm_BS_utilities(channel, network, new_partition, cell_assignment, I)

        # Complexity metrics
        state.num_utility_calculations += K
        state.num_longterm_rate_calculations += length(merged_block)
    end

    # Order the potential mergers
    ordered_mergers = Int[]
    if search_order == :greedy
        ordered_mergers = sortperm(squeeze(sum(merged_BS_utilities, 1), 1), rev=true)
    elseif search_order == :lexicographic
        ordered_mergers = 1:num_new_partitions
    end

    # Try to merge
    for sort_idx in ordered_mergers
        # Merge coalitions if everybody benefits
        if all(merged_BS_utilities[merged_BSs[sort_idx],sort_idx] .>= state.BS_utilities[merged_BSs[sort_idx]])
            state.partition = new_partitions[sort_idx]
            state.BS_utilities = merged_BS_utilities[:,sort_idx]
            state.r = min(length(new_partitions[sort_idx]), max_num_merging_coalitions)

            return true
        end
    end

    return false
end

##########################################################################
# BS utility definition
#
# Calculates the sum utility of the served MSs for each BS.
# cell_assignment and I are sent as part of the function signature to speed up
# evaluation slightly.
function longterm_BS_utilities(channel, network, partition, cell_assignment, I)
    BS_utilities = zeros(Float64, I)
    utilities, _ = longterm_utilities(channel, network, partition)
    for j = 1:I; for l in served_MS_ids(j, cell_assignment)
        BS_utilities[j] += utilities[l]
    end; end
    return BS_utilities
end
