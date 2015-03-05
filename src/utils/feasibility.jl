# Feasibility criterion from Liu2013
function is_IA_feasible(network, partition)
    Ns = get_no_MS_antennas(network)
    Ms = get_no_BS_antennas(network)
    ds = get_no_streams(network)
    assignment = get_assignment(network)

    # Can we even use these results?
    all(ds .== ds[1]) || Lumberjack.error("Feasibility check with Liu2013 only handles cases with equal number of streams per user.")
    d = ds[1]
    all(mod(Ns, d) .== 0) || Lumberjack.error("Feasibility check with Liu2013 only handles cases where d | N_k for all k.")
    all(mod(Ms, d) .== 0) || Lumberjack.error("Feasibility check with Liu2013 only handles cases where d | M_i for all i.")

    # List of served MSs
    served = [ served_MS_ids(i, assignment) for i = 1:length(Ms) ]
    served_length = [ length(served[i]) for i = 1:length(Ms) ]

    if all(Ns .== Ns[1]) && all(Ms .== Ms[1]) && all(served_length .== served_length[1])
        return Liu2013_IBC_symmetric(partition, Ns[1], Ms[1], d, served_length[1])
    else
        return Liu2013_IBC_heterogeneous(partition, Ns, Ms, d, served, served_length)
    end
end

function Liu2013_IBC_symmetric(partition, N, M, d, Kc)
    for block in partition.blocks
        I = length(block.elements)

        return ((I-1)*Kc*d <= M - Kc*d + N - d ? true : false)
    end
end

# Warning: This function is very slow.
function Liu2013_IBC_heterogeneous(partition, Ns, Ms, d, served, served_length)
    for block in partition.blocks
        # Condition (14a)
        for j in block.elements
            for i in block.elements; for k in served[i]
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
