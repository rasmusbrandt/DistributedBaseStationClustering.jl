##########################################################################
# Base station clustering based on coalitional formation. Three solution
# concepts are used: swap, individual and group-based stability.

##########################################################################
# Swap-based stability algorithm
type CoalitionFormationClustering_SwapState
    partition::Partition
    BS_throughputs::Vector{Float64}
    history::Vector{Set{IntSet}}
    num_searches::Vector{Int}
    num_sum_throughput_calculations::Int
end

function CoalitionFormationClustering_Swap(channel, network)
    I = get_no_BSs(network); K = get_no_MSs(network)

    aux_params = get_aux_assignment_params(network)
    @defaultize_param! aux_params "CoalitionFormationClustering_Swap:search_budget" 100
    search_budget = aux_params["CoalitionFormationClustering_Swap:search_budget"]
    @defaultize_param! aux_params "CoalitionFormationClustering_Swap:search_order" :random
    search_order = aux_params["CoalitionFormationClustering_Swap:search_order"]
    @defaultize_param! aux_params "CoalitionFormationClustering_Swap:stability_type" :swapee
    stability_type = aux_params["CoalitionFormationClustering_Swap:stability_type"]
    @defaultize_param! aux_params "CoalitionFormationClustering_Swap:use_history" true
    use_history = aux_params["CoalitionFormationClustering_Swap:use_history"]
    @defaultize_param! aux_params "CoalitionFormationClustering_Swap:starting_point" :singletons
    starting_point = aux_params["CoalitionFormationClustering_Swap:starting_point"]

    # Perform cell selection
    LargeScaleFadingCellAssignment!(channel, network)
    temp_cell_assignment = get_assignment(network)

    # Local RNG
    local_rng = MersenneTwister(int(time()))

    # Initial coalition structure
    if starting_point == :grand
        initial_partition = Partition(zeros(I))
    elseif starting_point == :singletons
        initial_partition = Partition(collect(0:(I-1)))
    end
    initial_BS_throughputs = longterm_BS_throughputs(channel, network, initial_partition, temp_cell_assignment, I)
    initial_num_searches = zeros(Int, I)
    state = CoalitionFormationClustering_SwapState(initial_partition, initial_BS_throughputs, [ Set{IntSet}() for i = 1:I ], initial_num_searches, K)

    # Let each BS deviate, and stop when no BS deviates (swap-based stability)
    deviation_performed = trues(I) # temporary, to enter the loop
    while any(deviation_performed)
        # Give all BSs a chance to deviate. If search_order == :greedy, we
        # let the BSs deviate in the order of their current throughputs, i.e. the
        # BS doing the best is going first. If search_order == :fair instead,
        # the BS which is doing the worst will go first. This is more similar
        # to GreedyClustering, where the strongest interfering links are clustered
        # first.
        if search_order == :greedy
            ordered_BS_list = sortperm(state.BS_throughputs, rev=true)
        elseif search_order == :fair
            ordered_BS_list = sortperm(state.BS_throughputs, rev=false)
        elseif search_order == :lexicographic
            ordered_BS_list = 1:I
        elseif search_order == :random
            ordered_BS_list = randperm(I, local_rng)
        end

        deviation_performed = falses(I)
        for i in ordered_BS_list
            deviation_performed[i] = deviate!(state, i, I, K, search_budget,
                stability_type, use_history, channel, network, temp_cell_assignment)
            #readline()
        end
    end
    throughputs, _, _, prelogs = longterm_throughputs(channel, network, state.partition)
    a = restricted_growth_string(state.partition)
    Lumberjack.info("CoalitionFormationClustering_Swap finished.",
        { :sum_throughput => sum(throughputs),
          :num_searches => state.num_searches,
          :a => a }
    )

    # Store prelogs for precoding
    set_aux_network_param!(network, prelogs[1], "prelogs_cluster_sdma")
    set_aux_network_param!(network, prelogs[2], "prelogs_network_sdma")

    # Store cluster assignment together with existing cell assignment
    network.assignment = Assignment(temp_cell_assignment.cell_assignment, cluster_assignment_matrix(network, state.partition))

    # Return results
    results = AssignmentResults()
    results["throughputs"] = throughputs
    results["a"] = a
    results["num_clusters"] = 1 + maximum(a)
    results["avg_cluster_size"] = avg_cluster_size(a)
    results["num_sum_throughput_calculations"] = state.num_sum_throughput_calculations
    results["num_searches"] = state.num_searches
    return results
end

# Lets BS i deviate in the swap stability model.
# Returns true if it did deviate, otherwise false.
function deviate!(state::CoalitionFormationClustering_SwapState, i, I, K,
    search_budget, stability_type, use_history, channel, network, cell_assignment)

    # println("BS $i is looking to deviate.")

    # First check that we have not exceeded our search budget
    if state.num_searches[i] >= search_budget
        return false
    end

    # Divide blocks such that my_old_block is the block that BS i used to
    # belong to, and other_blocks is an array of all other blocks.
    my_old_block = Block() # to store variable from inside loop
    other_blocks = Block[]
    for block in state.partition.blocks
        if i in block.elements
            my_old_block = Block(setdiff(block.elements, IntSet(i)))
        else
            push!(other_blocks, block)
        end
    end
    BSs_in_my_old_block = collect(my_old_block.elements)

    # Was this BS in a singleton coalition before?
    BS_in_singleton_coalition_before = length(my_old_block) == 0 ? true : false

    # Create all possible deviations for BS i
    Nother = 0
    for swapee_block in other_blocks
        if length(swapee_block) == 1 && BS_in_singleton_coalition_before
            continue
        else
            Nother += length(swapee_block.elements)
        end
    end
    num_new_partitions = length(other_blocks) + Nother + int(!BS_in_singleton_coalition_before)
    new_partitions = Array(Partition, num_new_partitions)
    deviated_BS_throughputs = zeros(Float64, I, num_new_partitions)
    swapees = zeros(Int64, num_new_partitions)
    local n

    # Deviations where BS i joins an existing coalition
    for n = 1:length(other_blocks)
        new_partition = Partition()
        other_blocks_cp = deepcopy(other_blocks) # see below
        union!(new_partition.blocks, other_blocks_cp)

        # Add the old block unless it used to be a singleton
        if !BS_in_singleton_coalition_before
            push!(new_partition.blocks, my_old_block)
        end

        # Add BS i to coalition n (this is why deepcopy is needed)
        push!(other_blocks_cp[n].elements, i)
        new_partitions[n] = new_partition
        deviated_BS_throughputs[:,n] = longterm_BS_throughputs(channel, network, new_partition, cell_assignment, I)

        # Complexity metrics
        state.num_sum_throughput_calculations += 1
    end

    # Deviations where BS i swaps with somebody in a full coalition
    for swapee_block in other_blocks
        # Two singletons should never swap
        if length(swapee_block) == 1 && BS_in_singleton_coalition_before
            continue
        end

        outside_blocks = Block[]
        for other_block2 in other_blocks
            if other_block2 == swapee_block; continue; end
            push!(outside_blocks, other_block2)
        end
        for j in swapee_block.elements
            new_partition = Partition()
            union!(new_partition.blocks, outside_blocks)

            # Add the old block unless it used to be a singleton
            if !BS_in_singleton_coalition_before
                push!(new_partition.blocks, my_old_block)
            end

            # Add BS i to and remove BS j from new block
            push!(new_partition.blocks, Block(union(setdiff(swapee_block.elements, IntSet(j)), IntSet(i))))

            # Add BS j to new singleton
            push!(new_partition.blocks, Block(IntSet(j)))

            n += 1
            new_partitions[n] = new_partition
            deviated_BS_throughputs[:,n] = longterm_BS_throughputs(channel, network, new_partition, cell_assignment, I)
            swapees[n] = j

            # Complexity metrics
            state.num_sum_throughput_calculations += 1
        end
    end

    # If BS i was not in a non-singleton coalition before deviation, add the
    # possibility that it belongs to a non-singleton coalition after deviation.
    if !BS_in_singleton_coalition_before
        new_partition = Partition()
        union!(new_partition.blocks, other_blocks)

        # Add the old block
        push!(new_partition.blocks, my_old_block)

        # Add BS i to new singleton coalition
        push!(new_partition.blocks, Block(IntSet(i)))
        new_partitions[end] = new_partition
        deviated_BS_throughputs[:,end] = longterm_BS_throughputs(channel, network, new_partition, cell_assignment, I)

        # Complexity metrics
        state.num_sum_throughput_calculations += 1
    end

    # println("$num_new_partitions number of possible new partitions:")
    # for partition in new_partitions
    #     println(partition)
    # end
    # println("With associated utilities")
    # println(deviated_BS_throughputs)

    # Check deviations, trying to join the coalitions in the order that
    # benefits BS i the most.
    for sort_idx in sortperm(squeeze(deviated_BS_throughputs[i,:], 1), rev=true)
        # Find block that BS i belongs to in this partition
        my_new_block = Block()
        for block in new_partitions[sort_idx].blocks
            if i in block.elements
                my_new_block = block
                break
            end
        end
        BSs_in_my_new_block = collect(my_new_block.elements)

        # Stop searching if we otherwise would exceed our search budget.
        if state.num_searches[i] >= search_budget
            return false
        end

        # Check that this coalition does not exist in the history, unless it is singleton
        if use_history && length(BSs_in_my_new_block) > 1 && IntSet(BSs_in_my_new_block) in state.history[i]
            return false
        end

        # Stop searching if we do not gain anymore.
        if deviated_BS_throughputs[i,sort_idx] <= state.BS_throughputs[i]
            return false
        end

        # Let's try to deviate
        state.num_searches[i] += 1

        # println("Wanting to join $my_new_block, resulting in partition $(new_partitions[sort_idx])")
        # if swapees[sort_idx] > 0
        #     println("===> (This means that $(swapees[sort_idx]) is being swapped out.)")
        # end

        # Check stability criterion
        if swap_stability(deviated_BS_throughputs[:,sort_idx], state.BS_throughputs, i, BSs_in_my_new_block, BSs_in_my_old_block, swapees[sort_idx], stability_type)
            # Let BS i join this coalition
            state.partition = new_partitions[sort_idx]
            state.BS_throughputs = deviated_BS_throughputs[:,sort_idx]

            # Add coalition to history
            push!(state.history[i], IntSet(BSs_in_my_new_block))

            # println("Join accepted")
            return true
        end
    end
    return false
end

# Check stability of a particular deviating BS for the swap
# coalition formation algorithm.
function swap_stability(new_BS_throughputs, old_BS_throughputs,
    deviating_BS_idx, new_coalition_idxs, old_coalition_idxs, swapee_BS_idx, stability_type)

    # Check that BS i strictly improves
    nash = (new_BS_throughputs[deviating_BS_idx] > old_BS_throughputs[deviating_BS_idx])
    if stability_type == :nash
        return nash
    end

    # Check if the BSs in the new coalition improve
    individual = (nash && all(new_BS_throughputs[new_coalition_idxs] .>= old_BS_throughputs[new_coalition_idxs]))
    if stability_type == :individual
        return individual
    end

    # Check if the swapee improves
    if swapee_BS_idx == 0
        # No swap happened
        swapee = individual
    else
        swapee = (individual && new_BS_throughputs[swapee_BS_idx] >= old_BS_throughputs[swapee_BS_idx])
    end
    if stability_type == :swapee
        return swapee
    end

    # Check if the BSs in the old coalition improve
    contractual = (swapee && all(new_BS_throughputs[old_coalition_idxs] .>= old_BS_throughputs[old_coalition_idxs]))
    if stability_type == :contractual
        return contractual
    end
end

##########################################################################
# Individual-based stability algorithm
type CoalitionFormationClustering_IndividualState
    partition::Partition
    BS_throughputs::Vector{Float64}
    history::Vector{Set{IntSet}}
    num_searches::Vector{Int}
    num_sum_throughput_calculations::Int
end

function CoalitionFormationClustering_Individual(channel, network)
    I = get_no_BSs(network); K = get_no_MSs(network)

    aux_params = get_aux_assignment_params(network)
    @defaultize_param! aux_params "CoalitionFormationClustering_Individual:search_budget" 100
    search_budget = aux_params["CoalitionFormationClustering_Individual:search_budget"]
    @defaultize_param! aux_params "CoalitionFormationClustering_Individual:search_order" :random
    search_order = aux_params["CoalitionFormationClustering_Individual:search_order"]
    @defaultize_param! aux_params "CoalitionFormationClustering_Individual:stability_type" :individual
    stability_type = aux_params["CoalitionFormationClustering_Individual:stability_type"]
    @defaultize_param! aux_params "CoalitionFormationClustering_Individual:use_history" true
    use_history = aux_params["CoalitionFormationClustering_Individual:use_history"]
    @defaultize_param! aux_params "CoalitionFormationClustering_Individual:starting_point" :singletons
    starting_point = aux_params["CoalitionFormationClustering_Individual:starting_point"]

    # Perform cell selection
    LargeScaleFadingCellAssignment!(channel, network)
    temp_cell_assignment = get_assignment(network)

    # Local RNG
    local_rng = MersenneTwister(int(time()))

    # Initial coalition structure
    if starting_point == :grand
        initial_partition = Partition(zeros(I))
    elseif starting_point == :singletons
        initial_partition = Partition(collect(0:(I-1)))
    end
    initial_BS_throughputs = longterm_BS_throughputs(channel, network, initial_partition, temp_cell_assignment, I)
    initial_num_searches = zeros(Int, I)
    state = CoalitionFormationClustering_IndividualState(initial_partition, initial_BS_throughputs, [ Set{IntSet}() for i = 1:I ], initial_num_searches, K)

    # Let each BS deviate, and stop when no BS deviates (individual-based stability)
    deviation_performed = trues(I) # temporary, to enter the loop
    while any(deviation_performed)
        # Give all BSs a chance to deviate. If search_order == :greedy, we
        # let the BSs deviate in the order of their current throughputs, i.e. the
        # BS doing the best is going first. If search_order == :fair instead,
        # the BS which is doing the worst will go first. This is more similar
        # to GreedyClustering, where the strongest interfering links are clustered
        # first.
        if search_order == :greedy
            ordered_BS_list = sortperm(state.BS_throughputs, rev=true)
        elseif search_order == :fair
            ordered_BS_list = sortperm(state.BS_throughputs, rev=false)
        elseif search_order == :lexicographic
            ordered_BS_list = 1:I
        elseif search_order == :random
            ordered_BS_list = randperm(I, local_rng)
        end

        deviation_performed = falses(I)
        for i in ordered_BS_list
            deviation_performed[i] = deviate!(state, i, I, K, search_budget,
                stability_type, use_history, channel, network, temp_cell_assignment)
        end
    end
    throughputs, _, _, prelogs = longterm_throughputs(channel, network, state.partition)
    a = restricted_growth_string(state.partition)
    Lumberjack.info("CoalitionFormationClustering_Individual finished.",
        { :sum_throughput => sum(throughputs),
          :num_searches => state.num_searches,
          :a => a }
    )

    # Store prelogs for precoding
    set_aux_network_param!(network, prelogs[1], "prelogs_cluster_sdma")
    set_aux_network_param!(network, prelogs[2], "prelogs_network_sdma")

    # Store cluster assignment together with existing cell assignment
    network.assignment = Assignment(temp_cell_assignment.cell_assignment, cluster_assignment_matrix(network, state.partition))

    # Return results
    results = AssignmentResults()
    results["throughputs"] = throughputs
    results["a"] = a
    results["num_clusters"] = 1 + maximum(a)
    results["avg_cluster_size"] = avg_cluster_size(a)
    results["num_sum_throughput_calculations"] = state.num_sum_throughput_calculations
    results["num_searches"] = state.num_searches
    return results
end

# Lets BS i deviate in the individual stability model.
# Returns true if it did deviate, otherwise false.
function deviate!(state::CoalitionFormationClustering_IndividualState, i, I, K,
    search_budget, stability_type, use_history, channel, network, cell_assignment)

    # First check that we have not exceeded our search budget
    if state.num_searches[i] >= search_budget
        return false
    end

    # Divide blocks such that my_old_block is the block that BS i used to
    # belong to, and other_blocks is an array of all other blocks.
    my_old_block = Block() # to store variable from inside loop
    other_blocks = Block[]
    for block in state.partition.blocks
        if i in block.elements
            my_old_block = Block(setdiff(block.elements, IntSet(i))) # BS i does not belong to the my_old_block anymore
        else
            push!(other_blocks, block)
        end
    end
    BSs_in_my_old_block = collect(my_old_block.elements)

    # Create all possible deviations for BS i
    BS_in_singleton_coalition_before = length(my_old_block) == 0 ? true : false # was this BS not in a singleton coalition before?
    num_new_partitions = length(other_blocks) + int(!BS_in_singleton_coalition_before)
    new_partitions = Array(Partition, num_new_partitions)
    deviated_BS_throughputs = zeros(Float64, I, num_new_partitions)

    # Deviations where BS i joins an existing coalition
    for n = 1:length(other_blocks)
        new_partition = Partition()
        other_blocks_cp = deepcopy(other_blocks) # see below
        union!(new_partition.blocks, other_blocks_cp)

        # Add the old block unless it used to be a singleton
        if !BS_in_singleton_coalition_before
            push!(new_partition.blocks, my_old_block)
        end

        # Add BS i to coalition n (this is why deepcopy is needed)
        push!(other_blocks_cp[n].elements, i)
        new_partitions[n] = new_partition
        deviated_BS_throughputs[:,n] = longterm_BS_throughputs(channel, network, new_partition, cell_assignment, I)

        # Complexity metrics
        state.num_sum_throughput_calculations += 1
    end

    # If BS i was in a non-singleton coalition before deviation, add the
    # possibility that it belongs to a non-singleton coalition after deviation.
    if !BS_in_singleton_coalition_before
        new_partition = Partition()
        union!(new_partition.blocks, other_blocks)

        # Add the old block
        push!(new_partition.blocks, my_old_block)

        # Add BS i to new singleton coalition
        push!(new_partition.blocks, Block(IntSet(i)))
        new_partitions[end] = new_partition
        deviated_BS_throughputs[:,end] = longterm_BS_throughputs(channel, network, new_partition, cell_assignment, I)

        # Complexity metrics
        state.num_sum_throughput_calculations += 1
    end

    # Check deviations, trying to join the coalitions in the order that
    # benefits BS i the most.
    for sort_idx in sortperm(squeeze(deviated_BS_throughputs[i,:], 1), rev=true)
        # Find block that BS i belongs to in this partition
        my_new_block = Block()
        for block in new_partitions[sort_idx].blocks
            if i in block.elements
                my_new_block = block
                break
            end
        end
        BSs_in_my_new_block = collect(my_new_block.elements)

        # Stop searching if we otherwise would exceed our search budget.
        if state.num_searches[i] >= search_budget
            return false
        end

        # Check that this coalition does not exist in the history, unless it is singleton
        if use_history && length(BSs_in_my_new_block) > 1 && IntSet(BSs_in_my_new_block) in state.history[i]
            return false
        end

        # Stop searching if we do not gain anymore.
        if deviated_BS_throughputs[i,sort_idx] <= state.BS_throughputs[i]
            return false
        end

        # Let's try to deviate
        state.num_searches[i] += 1

        # Check stability criterion
        if individual_stability(deviated_BS_throughputs[:,sort_idx], state.BS_throughputs, i, BSs_in_my_new_block, BSs_in_my_old_block, stability_type)
            # Let BS i join this coalition
            state.partition = new_partitions[sort_idx]
            state.BS_throughputs = deviated_BS_throughputs[:,sort_idx]

            # Add coalition to history
            push!(state.history[i], IntSet(BSs_in_my_new_block))

            return true
        end
    end
    return false
end

# Check stability of a particular deviating BS for the individual
# coalition formation algorithm.
function individual_stability(new_BS_throughputs, old_BS_throughputs,
    deviating_BS_idx, new_coalition_idxs, old_coalition_idxs, stability_type)

    # Check that BS i strictly improves
    nash = (new_BS_throughputs[deviating_BS_idx] > old_BS_throughputs[deviating_BS_idx])
    if stability_type == :nash
        return nash
    end

    # Check if the BSs in the new coalition improve
    individual = (nash && all(new_BS_throughputs[new_coalition_idxs] .>= old_BS_throughputs[new_coalition_idxs]))
    if stability_type == :individual
        return individual
    end

    # Check if the BSs in the old coalition improve
    contractual = (individual && all(new_BS_throughputs[old_coalition_idxs] .>= old_BS_throughputs[old_coalition_idxs]))
    if stability_type == :contractual
        return contractual
    end
end

##########################################################################
# Group-based stability algorithm
type CoalitionFormationClustering_GroupState
    partition::Partition
    BS_throughputs::Vector{Float64}
    r::Int
    num_sum_throughput_calculations::Int
end

function CoalitionFormationClustering_Group(channel, network)
    I = get_no_BSs(network); K = get_no_MSs(network)

    aux_params = get_aux_assignment_params(network)
    @defaultize_param! aux_params "CoalitionFormationClustering_Group:max_num_merging_coalitions" 3
    max_num_merging_coalitions = aux_params["CoalitionFormationClustering_Group:max_num_merging_coalitions"]
    @defaultize_param! aux_params "CoalitionFormationClustering_Group:search_order" :greedy
    search_order = aux_params["CoalitionFormationClustering_Group:search_order"]
    in(search_order, [:greedy, :lexicographic]) || Lumberjack.error("Incorrect CoalitionFormationClustering_Group:search_order.")

    # Perform cell selection
    LargeScaleFadingCellAssignment!(channel, network)
    temp_cell_assignment = get_assignment(network)

    # Initial coalition structure is the non-cooperative state
    initial_partition = Partition(collect(0:(I-1)))
    initial_BS_throughputs = longterm_BS_throughputs(channel, network, initial_partition, temp_cell_assignment, I)
    state = CoalitionFormationClustering_GroupState(initial_partition, initial_BS_throughputs, min(I, aux_params["CoalitionFormationClustering_Group:max_num_merging_coalitions"]), K)

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
    throughputs, _, _, prelogs = longterm_throughputs(channel, network, state.partition)
    a = restricted_growth_string(state.partition)
    Lumberjack.info("CoalitionFormationClustering_Group finished.",
        { :sum_throughput => sum(throughputs),
          :a => a }
    )

    # Store prelogs for precoding
    set_aux_network_param!(network, prelogs[1], "prelogs_cluster_sdma")
    set_aux_network_param!(network, prelogs[2], "prelogs_network_sdma")

    # Store cluster assignment together with existing cell assignment
    network.assignment = Assignment(temp_cell_assignment.cell_assignment, cluster_assignment_matrix(network, state.partition))

    # Return results
    results = AssignmentResults()
    results["throughputs"] = throughputs
    results["a"] = a
    results["num_clusters"] = 1 + maximum(a)
    results["avg_cluster_size"] = avg_cluster_size(a)
    results["num_sum_throughput_calculations"] = state.num_sum_throughput_calculations
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
    merged_BS_throughputs = zeros(Float64, I, num_new_partitions)
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
        merged_BS_throughputs[:,n] = longterm_BS_throughputs(channel, network, new_partition, cell_assignment, I)

        # Complexity metrics
        state.num_sum_throughput_calculations += 1
    end

    # Order the potential mergers
    ordered_mergers = Int[]
    if search_order == :greedy
        ordered_mergers = sortperm(squeeze(sum(merged_BS_throughputs, 1), 1), rev=true)
    elseif search_order == :lexicographic
        ordered_mergers = 1:num_new_partitions
    end

    # Try to merge
    for sort_idx in ordered_mergers
        # Merge coalitions if everybody benefits
        if all(merged_BS_throughputs[merged_BSs[sort_idx],sort_idx] .>= state.BS_throughputs[merged_BSs[sort_idx]])
            state.partition = new_partitions[sort_idx]
            state.BS_throughputs = merged_BS_throughputs[:,sort_idx]
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
function longterm_BS_throughputs(channel, network, partition, cell_assignment, I)
    BS_throughputs = zeros(Float64, I)
    throughputs, = longterm_throughputs(channel, network, partition)
    for j = 1:I; for l in served_MS_ids(j, cell_assignment)
        BS_throughputs[j] += throughputs[l]
    end; end
    return BS_throughputs
end

function randperm(n::Integer, rng)
    a = Array(typeof(n), n)
    a[1] = 1
    for i = 2:n
        j = iceil(i*rand(rng))
        a[i] = a[j]
        a[j] = i
    end
    return a
end
