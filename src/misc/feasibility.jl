##########################################################################
# IA feasibility check functions.
#
# These are based off of
# Liu, Yang, "On the Feasibility of Linear Interference Alignment for MIMO
# Interference Broadcast Channels With Constant Coefficients", IEEE Trans
# Signal Process, vol. 61, no. 9, pp. 2178-2191, 2013.
#

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
    # Check if we can use these results
    check_Liu2013_applicability(network)

    # Check IA feasibility with the tighest condition
    if check_Liu2013_symmetry(network)
        return Liu2013_IBC_symmetric(network, block)
    else
        return Liu2013_IBC_heterogeneous(network, block)
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
        Ns = get_num_MS_antennas(network)
        Ms = get_num_BS_antennas(network)
        ds = get_num_streams(network)

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

# Check if we can use the symmetric version of the condition
function check_Liu2013_symmetry(network)
    # Only check symmetry if we have not already done so
    if !haskey(network.aux_network_params, "check_Liu2013_symmetry:symmetry")
        I = get_num_BSs(network)
        Ns = get_num_MS_antennas(network)
        Ms = get_num_BS_antennas(network)
        ds = get_num_streams(network) # we already know that all MSs are served the same number of streams
        assignment = get_assignment(network)

        # Number of MSs served by each BS
        served_length = [ length(served_MS_ids(i, assignment)) for i = 1:I ]

        # Check symmetry and store in network, for caching purposes
        if all(Ns .== Ns[1]) && all(Ms .== Ms[1]) && all(served_length .== served_length[1])
            set_aux_network_param!(network, true, "check_Liu2013_symmetry:symmetry")

            # We also store all parameters needed to check feasibility, since
            # we are performing a huge number of feasibility checks in e.g.
            # exhaustive search or branch and bound.
            set_aux_network_param!(network, Ns[1], "check_Liu2013_symmetry:N")
            set_aux_network_param!(network, Ms[1], "check_Liu2013_symmetry:M")
            set_aux_network_param!(network, ds[1], "check_Liu2013_symmetry:d")
            set_aux_network_param!(network, served_length[1], "check_Liu2013_symmetry:Kc")
        else
            set_aux_network_param!(network, false, "check_Liu2013_symmetry:symmetry")
        end
    end

    return get_aux_network_param(network, "check_Liu2013_symmetry:symmetry")
end

# Special case for symmetric IBCs.
function Liu2013_IBC_symmetric(network, block)
    I = length(block)
    N = get_aux_network_param(network, "check_Liu2013_symmetry:N")
    M = get_aux_network_param(network, "check_Liu2013_symmetry:M")
    d = get_aux_network_param(network, "check_Liu2013_symmetry:d")
    Kc = get_aux_network_param(network, "check_Liu2013_symmetry:Kc")

    if I*Kc*d <= M + N - d
        return true
    else
        return false
    end
end

# General feasibility check for heterogeneous IBCs.
# Warning: This function is very slow since it checks all subsets of
# interfering links, which grows exponentially.
function Liu2013_IBC_heterogeneous(network, block)
    Lumberjack.warn("This function is very slow, since it checks all subsets of interfering links. Perhaps a polynomial check, similar to de Kerret for the IC, could be used?")

    I = get_num_BSs(network)
    Ns = get_num_MS_antennas(network); Ms = get_num_BS_antennas(network)
    d = get_num_streams(network)[1]
    served = [ served_MS_ids(i, assignment) for i = 1:I ]
    served_length = [ length(served[i]) for i = 1:I ]

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
