simulation_params = [
    "I" => 25, "Kc" => 1, "N" => 2, "M" => 4, "d" => 1,
    "Ndrops" => 10, "Nsim" => 1,
    "geography_size" => (250.,250.),
    "MS_serving_BS_distance" => 50.,
    "assignment_methods" => [
        # ExhaustiveSearchClustering,
        # BranchAndBoundClustering,

        CoalitionFormationClustering_Group,
        CoalitionFormationClustering_Individual,

        GreedyClustering_Single,
        GreedyClustering_Multiple,

        # Chen2014_ExhaustiveSearch,
        # Peters2012_Heuristic,

        GrandCoalitionClustering,
        RandomClustering,
        NoClustering,
    ],
    "aux_assignment_params" => [
        "max_MSs_per_BS" => 1,

        "clustering_type" => :spectrum_sharing,
        "apply_overhead_prelog" => true,
        "IA_infeasible_negative_inf_utility" => true,
        "replace_E1_utility_with_lower_bound" => false,

        "BranchAndBoundClustering:bracket_E1" => false,

        "CoalitionFormationClustering_Group:max_merge_size" => 3,
        "CoalitionFormationClustering_Group:search_order" => :greedy,
        "CoalitionFormationClustering_Individual:search_budget" => 10,
        "CoalitionFormationClustering_Individual:search_order" => :greedy,
        "CoalitionFormationClustering_Individual:stability_type" => :contractual,
    ],
    "aux_precoding_params" => [
        "initial_precoders" => "eigendirection",
        "stop_crit" => 1e-2,
        "max_iters" => 1000,
    ],
]
