function is_IA_feasible(partition, Ns, Ms, ds)
    all(ds .== ds[1]) || Lumberjack.error("IA feasibility only implemented for symmetric d case.")

    for block in partition.blocks
        # Properness check
        I = length(block)
        if I*(I+1)*ds[1] > sum(Ms[block.elements]) + sum(Ns[block.elements])
            return false
        end
    end
    return true
end
