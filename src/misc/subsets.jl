##########################################################################
# Helper functions for subset generation.

# Returns all non-empty subsets of supplied list of elements
function all_nonempty_subsets(elements::Vector)
    N = length(elements)

    subsets = Set[]

    for pattern = ((1 << N) - 1):-1:1
        subset = Set()
        for element_idx = 1:N
            # Add elements to subset based on pattern
            if ((pattern >> (element_idx - 1)) & 1) != 0
                push!(subset, elements[element_idx])
            end
        end
        push!(subsets, subset)
    end

    return subsets
end
