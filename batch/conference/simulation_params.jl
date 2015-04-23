fc = 5e9 # GHz
Wc = 300e3 # kHz

c = 300e6 # m/s
λ = c/fc # m

BS_density = 1/(1000*1000) # BSs per m^2

simulation_params = [
    "Ndrops" => 10, "Nsim" => 5,
    "MS_serving_BS_distance" => nothing, # random placement of MSs with greedy user association
    "aux_assignment_params" => [
        "max_MSs_per_BS" => 1,

        "clustering_type" => :spectrum_sharing,
        "apply_overhead_prelog" => true,
        "IA_infeasible_negative_inf_utility" => true,
        "replace_E1_utility_with_lower_bound" => false,

        "CoalitionFormationClustering_Individual:search_order" => :lexicographic,
        "CoalitionFormationClustering_Individual:stability_type" => :individual,
        "CoalitionFormationClustering_Individual:use_history" => true,
    ],
    "aux_precoding_params" => [
        "initial_precoders" => "eigendirection",
        "stop_crit" => 1e-3,
        "max_iters" => 1000,
    ],
]
