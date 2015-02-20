function feasibility_IA(partition, Ns, Ms, ds)
    for block in partition.blocks
        # Properness check
        K = length(block)
        if ds[1] > 1/(K*(K+1))*(sum(Ms) + sum(Ns))
            return false
        end

        # if length(block) > Ms[1] + Ns[1] - 1
        #     return false
        # end
    end
    return true
end
