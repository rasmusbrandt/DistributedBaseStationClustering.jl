function GrandCoalitionClustering(channel, network)
    K = get_no_MSs(network)

    assign_cells_by_id!(network)
    network.assignment = Assignment(network.assignment.cell_assignment, ones(K,K))
end
