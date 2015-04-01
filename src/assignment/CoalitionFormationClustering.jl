# Base station clustering based on coalitional formation.

# Individual-based stability.
function CoalitionFormationClustering_Individual(channel, network)
    I = get_no_BSs(network)

    aux_params = get_aux_assignment_params(network)
    @defaultize_param! aux_params "CoalitionFormationClustering:search_budget" 50
    @defaultize_param! aux_params "CoalitionFormationClustering:rng_seed" 7326253

    # Get our own rng for the randperm operation
    rng = MersenneTwister(aux_params["CoalitionFormationClustering:rng_seed"])

    # Perform cell selection
    LargeScaleFadingCellAssignment!(channel, network)
    temp_cell_assignment = get_assignment(network)

    # Initial coalition structure is the non-cooperative state
    partition = Partition(collect(0:(I-1)))
    current_BS_utilities = longterm_BS_utilities(channel, network, partition, temp_cell_assignment, I)

    # Search budgets
    no_searches = zeros(Int, I)

    # Perform coalitional formation until no BS deviates
    deviation_performed = true
    while deviation_performed
        deviation_performed = false

        # Give all BSs a chance to deviate. Do this in a random order.
        for i in randperm(rng, I)
            # Only continue with this BS if it has enough search budget left
            if no_searches[i] >= aux_params["CoalitionFormationClustering:search_budget"]
                continue # for i in randperm(rng, I)
            end

            # Divide blocks such that old_block is the block that BS i used to
            # belong to, and other_blocks is an array of all other blocks.
            old_block = Block() # to store variable from inside loop
            other_blocks = Block[]
            for block in partition.blocks
                if i in block.elements
                    old_block = Block(setdiff(block.elements, IntSet(i))) # BS i does not belong to the old_block anymore
                else
                    push!(other_blocks, block)
                end
            end

            # Create all possible deviations for BS i
            BS_grouped_before = length(old_block) > 0 ? true : false # was this BS not in a singleton coalition before?
            no_new_partitions = length(other_blocks) + int(BS_grouped_before)
            new_partitions = Array(Partition, no_new_partitions)
            deviating_BS_utilities = zeros(Float64, I, no_new_partitions)
            for n = 1:length(other_blocks)
                # Loop over the deviations where BS i joins an existing coalition
                new_partition = Partition()
                other_blocks_cp = deepcopy(other_blocks) # need to deepcopy, so the created coalitions will not all be the same...
                union!(new_partition.blocks, other_blocks_cp)

                # Add the old block unless it used to be a singleton
                if BS_grouped_before
                    push!(new_partition.blocks, old_block)
                end

                # Add BS i to coalition n
                push!(other_blocks_cp[n].elements, i)

                new_partitions[n] = new_partition
                deviating_BS_utilities[:,n] = longterm_BS_utilities(channel, network, new_partition, temp_cell_assignment, I)
            end
            if BS_grouped_before
                # BS i was in a non-singleton coalition before. Add the
                # the possibility that it belongs to a non-singleton
                # coalition after deviation.
                new_partition = Partition()
                other_blocks_cp = deepcopy(other_blocks)
                union!(new_partition.blocks, other_blocks_cp)

                # Add the old block
                push!(new_partition.blocks, old_block)

                # Add BS i to new singleton coalition
                push!(new_partition.blocks, Block(IntSet(i)))

                new_partitions[end] = new_partition
                deviating_BS_utilities[:,end] = longterm_BS_utilities(channel, network, new_partition, temp_cell_assignment, I)
            end

            # Check deviations, trying to join the coalitions in the order that
            # benefits BS i the most.
            for sort_idx in sortperm(squeeze(deviating_BS_utilities[i,:], 1), rev=true)
                # Do we even want to join this coalition?
                if deviating_BS_utilities[i,sort_idx] > current_BS_utilities[i]
                    # Do we have enough search budget left?
                    if no_searches[i] < aux_params["CoalitionFormationClustering:search_budget"]
                        no_searches[i] += 1

                        # Find block that BS i belongs to in this partition
                        my_block = Block()
                        for block in new_partitions[sort_idx].blocks
                            if i in block.elements
                                my_block = block
                            end
                        end

                        # Check if the existing members of this coalition allow BS i to join (this check includes BS i, unnecessarily)
                        BSs_in_block = collect(my_block.elements)
                        if all(deviating_BS_utilities[BSs_in_block,sort_idx] .> current_BS_utilities[BSs_in_block])
                            # Switch BS i to this coalition
                            partition = new_partitions[sort_idx]
                            current_BS_utilities = deviating_BS_utilities[:,sort_idx]

                            deviation_performed = true
                            break # for sort_idx in dev_BS_utils_sort
                        end
                    end
                else
                    # Not enough search budget left
                    break # for sort_idx in dev_BS_utils_sort
                end
            end
        end
    end
    utilities, _ = longterm_utilities(channel, network, partition)
    Lumberjack.info("CoalitionFormationClustering finished.", { :sum_utility => sum(utilities), :a => restricted_growth_string(partition), :no_searches => no_searches })

    # Store cluster assignment together with existing cell assignment
    temp_assignment = get_assignment(network)
    network.assignment = Assignment(temp_assignment.cell_assignment, cluster_assignment_matrix(network, partition))

    # Return results
    results = AssignmentResults()
    results["utilities"] = utilities
    return results
end

# Calculates the sum utility of the served MSs for each BS. assignment and I
# are sent as part of the function signature to speed up evaluation slightly.
function longterm_BS_utilities(channel, network, partition, assignment, I)
    BS_utilities = zeros(Float64, I)
    utilities, _ = longterm_utilities(channel, network, partition)
    for j = 1:I; for l in served_MS_ids(j, assignment)
        BS_utilities[j] += utilities[l]
    end; end
    return BS_utilities
end

# randperm taken from base, extended to accept a rng
function randperm(rng, n::Integer)
    a = Array(typeof(n), n)
    a[1] = 1
    for i = 2:n
        j = iceil(i*rand(rng))
        a[i] = a[j]
        a[j] = i
    end
    return a
end
