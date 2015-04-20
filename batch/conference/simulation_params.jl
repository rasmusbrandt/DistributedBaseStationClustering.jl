fc = 5e9 # GHz
Wc = 300e3 # kHz

c = 300e6 # m/s
λ = c/fc # m

simulation_params = [
    "Ndrops" => 100, "Nsim" => 5,
    "geography_size" => (250.,250.),
    "MS_serving_BS_distance" => nothing, # random placement of MSs with greedy user association
    "aux_assignment_params" => [
        "max_MSs_per_BS" => 1,

        "clustering_type" => :spectrum_sharing,
        "apply_overhead_prelog" => true,
        "IA_infeasible_negative_inf_utility" => true,
        "replace_E1_utility_with_lower_bound" => false,

        "CoalitionFormationClustering_Individual:search_budget" => 10,
        "CoalitionFormationClustering_Individual:search_order" => :greedy,
        "CoalitionFormationClustering_Individual:stability_type" => :contractual,
    ],
    "aux_precoding_params" => [
        "initial_precoders" => "dft",
        "stop_crit" => 1e-3,
        "max_iters" => 1000,
    ],
]
