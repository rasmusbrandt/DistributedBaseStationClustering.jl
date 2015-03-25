##########################################################################
# IA feasibility checks

# Returns true if IA is feasible for this partition
function is_IA_feasible(network, partition::Partition)
    check_Liu2013_applicability(network)

    for block in partition.blocks
        is_IA_feasible(network, block) || return false
    end
    return true
end

# Returns true if IA is feasible for this block
function is_IA_feasible(network, block::Block)
    check_Liu2013_applicability(network)

    I = get_no_BSs(network)
    Ns = get_no_MS_antennas(network)
    Ms = get_no_BS_antennas(network)
    ds = get_no_streams(network); d = ds[1]
    assignment = get_assignment(network)

    # List of served MSs
    served = [ served_MS_ids(i, assignment) for i = 1:I ]
    served_length = [ length(served[i]) for i = 1:I ]

    if all(Ns .== Ns[1]) && all(Ms .== Ms[1]) && all(served_length .== served_length[1])
        return Liu2013_IBC_symmetric(block, Ns[1], Ms[1], d, served_length[1])
    else
        return Liu2013_IBC_heterogeneous(block, Ns, Ms, d, served, served_length)
    end
end

# Check that the Liu2013 IA feasibility results can be applied
function check_Liu2013_applicability(network)
    if haskey(network.aux_network_params, "check_Liu2013_applicability:applicable")
        # If we have already checked the applicability, don't do it again.
        # If the results turned out to be not applicable, then an error has
        # been raised already.
        return true
    else
        Ns = get_no_MS_antennas(network)
        Ms = get_no_BS_antennas(network)
        ds = get_no_streams(network)

        # Can we even use these results?
        # This check is actually slightly more general than necessary, since
        # it checks the requirements for the entire network. In principle, we
        # could check feasibility for just this block.
        all(ds .== ds[1]) || Lumberjack.error("Feasibility check with Liu2013 only handles cases with equal number of streams per user.")
        d = ds[1]
        all(mod(Ns, d) .== 0) || Lumberjack.error("Feasibility check with Liu2013 only handles cases where d | N_k for all k.")
        all(mod(Ms, d) .== 0) || Lumberjack.error("Feasibility check with Liu2013 only handles cases where d | M_i for all i.")

        # Store in network, for caching purposes.
        set_aux_network_param!(network, true, "check_Liu2013_applicability:applicable")
    end
end

# Special case for symmetric IBCs.
function Liu2013_IBC_symmetric(block, N, M, d, Kc)
    I = length(block.elements)

    if I*Kc*d <= M + N - d
        return true
    else
        return false
    end
end

# General feasibility check for heterogeneous IBCs.
# Warning: This function is very slow since it checks all subsets of
# interfering links, which grows exponentially.
function Liu2013_IBC_heterogeneous(block, Ns, Ms, d, served, served_length)
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

    return true
end
