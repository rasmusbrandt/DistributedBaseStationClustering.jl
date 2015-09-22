simulation_params["assignment_methods"] = [ CoalitionFormationClustering_AttachOrSupplant, ]

simulation_params["precoding_methods"] = [
    RobustIntraclusterWMMSE,
    NaiveIntraclusterWMMSE,

    RobustChen2014_MaxSINR,
    NaiveChen2014_MaxSINR,
]
