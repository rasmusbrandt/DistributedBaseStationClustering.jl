function partition_to_assignment_matrix(partition, K)
    assignment_matrix = eye(K, K)

    # Loop over IA clusters, which are described by
    # the blocks in the set partition.
    for block in partition.blocks
        intercluster_matrix = zeros(Int, K, K)
        for k in block.elements
            for l in setdiff(block.elements, k)
                intercluster_matrix[k,l] = 1
            end
        end
        assignment_matrix += intercluster_matrix
    end

    return assignment_matrix
end
