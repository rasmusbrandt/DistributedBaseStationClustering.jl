simulation_params["assignment_methods"] = [ BranchAndBoundClustering, ]

simulation_params["precoding_methods"] = [
    RobustIntraclusterWMMSE,
    NaiveIntraclusterWMMSE,

    RobustChen2014_MaxSINR,
    NaiveChen2014_MaxSINR,
]
