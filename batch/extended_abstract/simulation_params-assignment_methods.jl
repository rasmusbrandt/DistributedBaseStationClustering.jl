simulation_params["assignment_methods"] = [
    # ExhaustiveSearchClustering,
    # BranchAndBoundClustering,

    # CoalitionFormationClustering_Group,
    CoalitionFormationClustering_Individual,

    # GreedyClustering_Single,
    GreedyClustering_Multiple,

    # Chen2014_ExhaustiveSearch,
    # Peters2012_Heuristic,

    GrandCoalitionClustering,
    RandomClustering,
    NoClustering,
]

simulation_params["precoding_methods"] = [ RobustIntraclusterWMMSE ]
