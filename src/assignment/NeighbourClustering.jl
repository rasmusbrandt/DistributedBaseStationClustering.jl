function NeighbourClustering(channel, network)
    assign_cells_by_id!(network)
    network.assignment = Assignment(network.assignment.cell_assignment,
            [1 0 0 1 1 1;
             0 1 0 1 1 1;
             0 0 1 1 1 1;
             1 1 1 1 0 0;
             1 1 1 0 1 0;
             1 1 1 0 0 1;]
    )
end
