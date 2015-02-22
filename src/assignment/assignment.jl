function partition_to_assignment_matrix(partition, K, I, assignment)
    assignment_matrix = zeros(Int, K, I)

    for block in partition.blocks
        for i in block.elements
            for j in block.elements; for l in served_MS_ids(j, assignment)
                assignment_matrix[l,i] = 1
            end; end
        end
    end

    return assignment_matrix
end
