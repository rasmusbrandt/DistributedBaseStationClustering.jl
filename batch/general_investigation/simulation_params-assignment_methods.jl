simulation_params["assignment_methods"] = [
    # ExhaustiveSearchClustering,
    # BranchAndBoundClustering,

    CoalitionFormationClustering_Swap,
    CoalitionFormationClustering_Individual,
    # CoalitionFormationClustering_Group,

    GreedyClustering_Single,
    GreedyClustering_Multiple,

    # Chen2014_LinearObj_ExhaustiveSearch,
    # Peters2012_Heuristic,

    GrandCoalitionClustering,
    RandomClustering,
    NoClustering,
]

simulation_params["precoding_methods"] = [ RobustIntraclusterWMMSE, ]
