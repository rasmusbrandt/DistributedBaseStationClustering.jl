##########################################################################
# Partitions

# A block is just a set of integers.
immutable Block
    elements::Vector{Int}
end
Block() = Block([])
getindex(b::Block, idx) = b.elements[idx]
Base.push!(b::Block, e::Int) = push!(b.elements, e)
Base.show(io::IO, b::Block) = showcompact(io, b.elements)
Base.length(b::Block) = length(b.elements)

# The Partition type describes a partition of integers into blocks.
immutable Partition
    blocks::Vector{Block}
end
Base.show(io::IO, p::Partition) = showcompact(io, p.blocks)
Base.length(p::Partition) = length(p.blocks)

# Create partition from restricted growth string
function Partition(rgs::Vector)
    no_blocks = maximum(rgs) + 1
    blocks = [ Block() for n = 1:no_blocks ]

    for l = 1:length(rgs)
        push!(blocks[rgs[l] + 1], l)
    end

    return Partition(blocks)
end

# PartitionIterator gives all partitions of the set 1:n. This iterator is
# based on Algorithm H from TAoCP 7.2.1.5.
immutable PartitionIterator
    n::Int

    PartitionIterator(n) = n >= 2 ? new(n) : error("only applicable for n >= 2.")
end

type PartitionIteratorState
    a::Vector{Int}
    b::Vector{Int}
    m::Int
    terminate::Bool
end
PartitionIteratorState(n) = PartitionIteratorState(zeros(n), ones(n - 1), 1, false)

Base.start(p::PartitionIterator) = PartitionIteratorState(p.n) # H1
Base.done(p::PartitionIterator, s) = s.terminate
function Base.next(p::PartitionIterator, s)
    # H2
    current_partition = Partition(s.a)

    if s.a[p.n] < s.m
        # H3
        s.a[p.n] += 1
    else
        # H4
        j = p.n - 1
        while s.a[j] == s.b[j]
            j -= 1
        end

        # H5
        if j == 1
            s.terminate = true
        else
            s.a[j] += 1

            # H6
            s.m = s.b[j] + (s.a[j] == s.b[j] ? 1 : 0)
            j += 1
            while j < p.n
                s.a[j] = 0
                s.b[j] = s.m
                j += 1
            end
            s.a[p.n] = 0
        end
    end

    # H2 cont'd.
    return current_partition, s
end

##########################################################################
# Stirling and Bell numbers

# Stirling2NumberCoefficientIterator gives all coefficients that when summed
# result in a Stirling number of the second kind. Note that the coefficients
# may be real numbers, but when summed are integers. The coefficients are
# thus calculated as floating point numbers.
immutable Stirling2NumberCoefficientIterator
    n::Int
    k::Int
    k_fac::Int
end
Stirling2NumberCoefficientIterator(n, k) = Stirling2NumberCoefficientIterator(n, k, factorial(k))
Base.start(s::Stirling2NumberCoefficientIterator) = 0 # state: k
Base.done(s::Stirling2NumberCoefficientIterator, j) = (j > s.k)
Base.next(s::Stirling2NumberCoefficientIterator, j) = ((iseven(s.k - j) ? 1. : -1.)*(binomial(s.k, j)*j^s.n/s.k_fac), j + 1)

# Stirling2NumberIterator gives all Stirling numbers of the second kind
# up to n. The Stirling numbers are integers, but the Stirling number
# coefficients are not. We therefore round the sum of the coefficients,
# in order to end up with an integer.
immutable Stirling2NumberIterator
    n::Int
end
Base.start(s::Stirling2NumberIterator) = 0 # state: k
Base.done(s::Stirling2NumberIterator, k) = (k > s.n)
Base.next(s::Stirling2NumberIterator, k) = (iround(sum(Stirling2NumberCoefficientIterator(s.n, k))), k + 1)

# The Bell number can we calculated as the sum of all
# Stirling numbers of the second kind.
bell_number(n) = sum(Stirling2NumberIterator(n))
