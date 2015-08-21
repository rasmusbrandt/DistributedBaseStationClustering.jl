##########################################################################
# Group-based stability algorithm
type CoalitionFormationClustering_GroupState
    partition::Partition
    BS_throughputs::Vector{Float64}
    r::Int
    num_sum_throughput_calculations::Int
end

function CoalitionFormationClustering_Group(channel, network)
    I = get_num_BSs(network); K = get_num_MSs(network)

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