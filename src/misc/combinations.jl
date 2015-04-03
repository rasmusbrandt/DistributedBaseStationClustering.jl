##########################################################################
# Helper functions to generate combinations.

# kCombinationIterator gives all k-combinations of the set 1:n. The iterator
# is based on
#
# Ruskey, "Combinatorial Generation", University of Victoria, 2003.
#
# Note that a is one-based here, contrarily to partitions.jl.
immutable kCombinationIterator
    n::Int
    k::Int

    kCombinationIterator(n, k) = (0 <= k <= n) ? new(n, k) : error("only applicable for 0 <= k <= n.")
end

type kCombinationIteratorState
    a::Vector{Int}
    terminate::Bool
end
kCombinationIteratorState(k) = kCombinationIteratorState(1:k, false)

Base.start(c::kCombinationIterator) = kCombinationIteratorState(c.k)
Base.done(c::kCombinationIterator, s) = s.terminate
function Base.next(c::kCombinationIterator, s)
    current_combination = copy(s.a)

    # C5
    j = c.k

    # C6
    while j > 0 && s.a[j] == c.n - c.k + j
        j -= 1
    end

    # C7
    if j == 0
        s.terminate = true
    else
        # C8
        s.a[j] += 1

        # C9
        for i = (j+1):c.k
            s.a[i] = s.a[i-1] + 1
        end
    end

    return current_combination, s
end
