# Implementation of feasibility criterion from Razaviyayn2012 and Liu2013
function is_IA_feasible(partition, Ns, Ms, ds, assignment)
    # Can we even use these results?
    all(ds .== ds[1]) || Lumberjack.error("Feasibility check only handles cases with equal number of streams per user.")
    d = ds[1]
    all(mod(Ns, d) .== 0) || Lumberjack.error("Feasibility check only handles cases where d | N_k for all k.")
    all(mod(Ms, d) .== 0) || Lumberjack.error("Feasibility check only handles cases where d | M_i for all i.")

    # List of served MSs
    served = [ served_MS_ids(i, assignment) for i = 1:length(Ms) ]
    served_length = [ length(served[i]) for i = 1:length(Ms) ]

    # IC or IBC?
    if all(served_length .<= 1)
        # Is this a symmetric scenario?
        if all(Ns .== Ns[1]) && all(Ms .== Ms[1])
            return Razaviyayn2012_IC_symmetric(partition, Ns[1], Ms[1], d)
        else
            # Implement this.
        end
    else
        # Is this a symmetric scenario?
        if all(Ns .== Ns[1]) && all(Ms .== Ms[1]) && all(served_length .== served_length[1])
            return Liu2013_IBC_symmetric(partition, Ns[1], Ms[1], d, served_length[1])
        else
            return Liu2013_IBC_heterogeneous(partition, Ns, Ms, d, assignment)
        end
    end
end

# IA feasibility for a symmetric IC
function Razaviyayn2012_IC_symmetric(partition, N, M, d)
    for block in partition.blocks
        K = length(block.elements)

        return (K*d <= M + N - 1 ? true : false)
    end
end

# IA feasibility for a symmetric IBC
function Liu2013_IC_symmetric(partition, N, M, d, Kc)
    for block in partition.blocks
        I = length(block.elements)

        return ((I-1)*Kc*d <= M - Kc*d + N - d ? true : false)
    end
end

# IA feasibility for a heterogeneous IBC
function Liu2013_IBC_heterogeneous(partition, Ns, Ms, d, assignment)
    for block in partition.blocks
        # Condition (14a)
        for j in block.elements
            for i in block.elements; for k in served_MS_ids(i, assignment)
                min(Ms[j] - served_length[j]*d, Ns[k] - d) >= 0 || return false
            end; end
        end

        # Enumerate all interfering links
        interfering_cells = (Int, Int)[]
        for i in block.elements; for j in block.elements
            if i != j
                push!(interfering_cells, (i, j))
            end
        end; end

        # Condition (14b)
        for interferering_cells_subset in all_nonempty_subsets(interfering_cells)
            is = unique([ i for (i, _) in interferering_cells_subset ])
            js = unique([ j for (_, j) in interferering_cells_subset ])

            term1 = 0
            for j in js
                term1 += (Ms[j] - served_length[j]*d)*served_length[j]
            end

            term2 = 0
            for i in is; for k in served[i]
                term2 += Ns[k] - d
            end; end

            term3 = 0
            for (i, j) in interferering_cells_subset
                term3 += served_length[j]*served_length[i]*d
            end
            term1 + term2 >= term3 || return false
        end
    end

    return true
end
