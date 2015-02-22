immutable Block
    elements::Vector{Int}
end
Block() = Block([])
Block(e::Int) = Block([e])
Block(e1::Int, e2::Int) = Block([e1, e2])
Block(elements::Vector{Int}, extra_element::Int) = Block(vcat(elements, extra_element))
getindex(b::Block, idx) = b.elements[idx]
Base.show(io::IO, b::Block) = showcompact(io, b.elements)
Base.push!(b::Block, e::Int) = push!(b.elements, e)
Base.length(b::Block) = length(b.elements)

immutable Partition
    blocks::Vector{Block}
end
Partition() = Partition([])
Partition(b::Block) = Partition([b])
Partition(b1::Block, b2::Block) = Partition([b1, b2])
Partition(blocks::Vector{Block}, extra_block::Block) = Partition(vcat(blocks, extra_block))
getindex(p::Partition, idx) = p.blocks[idx]
Base.show(io::IO, p::Partition) = showcompact(io, p.blocks)
Base.length(p::Partition) = length(p.blocks)

add_element_to_block!(p::Partition, e::Int, block_idx::Int) =
    push!(p[block_idx], e)

function all_twoblock_partitions(elements::Block)
    N = length(elements)

    # No result if less than two elements
    if N < 2
        return Partition[]
    end

    # Generate all two-block partitions
    partitions = Partition[]
    for pattern = 1:((1 << (N-1)) - 1)
        # Create new partition. We specifically put the first element in the
        # first subset, to avoid ordering problems.
        partition = Partition(Block(elements[1]), Block())

        # Distribute the remaining elements into the two blocks
        for idx = 2:N
            add_element_to_block!(partition, elements[idx], ((pattern >> (idx - 2)) & 1) + 1)
        end

        # Store this partition in the list of all partitions
        push!(partitions, partition)
    end

    return partitions
end

function all_partitions(elements)
    partitions = Partition[]
    all_partitions(partitions, Block[], Block(elements))
    return partitions
end

function all_partitions(partitions::Vector{Partition}, fixed_blocks::Vector{Block}, suffix_elements::Block)
    # The trivial partition is adding the suffix elements as one block
    # after all fixed blocks.
    push!(partitions, Partition(fixed_blocks, suffix_elements))

    # Get all two-block partitions of the suffix elements, and sub-divide
    # them recursively.
    suffix_partitions = all_twoblock_partitions(suffix_elements)
    sub_partitions_return = Partition[]
    for suffix_partition in suffix_partitions
        sub_partitions = all_partitions(partitions, push!(copy(fixed_blocks), suffix_partition[1]), suffix_partition[2])
        for sub_partition in sub_partitions
            push!(partitions, sub_partition)
            push!(sub_partitions_return, sub_partition)
        end
    end

    return sub_partitions_return
end
