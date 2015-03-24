# Creates a K-by-I matrix whose (i,k) entry denotes whether or not
# BS i coordinates to MS k. Note that this matrix is different than
# the partition logical matrix in partitions.jl, which is an
# I-by-I matrix.
function cluster_assignment_matrix(network, partition)
    I = get_no_BSs(network); K = get_no_MSs(network)
    assignment = get_assignment(network)

    A = zeros(Int, K, I)

    for block in partition.blocks
        for i in block.elements
            for j in block.elements; for l in served_MS_ids(j, assignment)
                A[l,i] = 1
            end; end
        end
    end

    return A
end

include("BranchAndBoundClustering.jl")
include("Chen2014_LinearObjClustering.jl")
include("ExhaustiveSearchClustering.jl")
include("GrandCoalitionClustering.jl")
include("GreedyClustering.jl")
include("NoClustering.jl")
include("RandomClustering.jl")
