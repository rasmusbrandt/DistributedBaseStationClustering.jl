function partition_to_cluster_assignment_matrix(partition, K, I, assignment)
    cluster_assignment_matrix = zeros(Int, K, I)

    for block in partition.blocks
        for i in block.elements
            for j in block.elements; for l in served_MS_ids(j, assignment)
                cluster_assignment_matrix[l,i] = 1
            end; end
        end
    end

    return cluster_assignment_matrix
end

include("Chen2014_LinearObjClustering.jl")
include("ExhaustiveSearchClustering.jl")
include("GrandCoalitionClustering.jl")
include("NoClustering.jl")
include("RandomClustering.jl")
