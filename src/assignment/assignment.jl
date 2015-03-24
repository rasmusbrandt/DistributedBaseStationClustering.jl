function cluster_assignment_matrix(network, partition)
    I = get_no_BSs(network); K = get_no_MSs(network)
    assignment = get_assignment(network)

    cl_assignment_matrix = zeros(Int, K, I)

    for block in partition.blocks
        for i in block.elements
            for j in block.elements; for l in served_MS_ids(j, assignment)
                cl_assignment_matrix[l,i] = 1
            end; end
        end
    end

    warn("Remove me.")

    return cl_assignment_matrix
end

include("BranchAndBoundClustering.jl")
include("Chen2014_LinearObjClustering.jl")
include("ExhaustiveSearchClustering.jl")
include("GrandCoalitionClustering.jl")
include("GreedyClustering.jl")
include("NoClustering.jl")
include("RandomClustering.jl")
