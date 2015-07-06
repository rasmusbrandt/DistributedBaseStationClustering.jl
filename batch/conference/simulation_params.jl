fc = 2e9 # GHz
Wc = 300e3 # kHz

c = 300e6 # m/s
λ = c/fc # m

ISD = 500
BS_density = 1/(sqrt(3)/2*ISD^2) # BSs per m^2, for hexagonal cells

simulation_params = [
    "Ndrops" => 100, "Nsim" => 5,
    "MS_serving_BS_distance" => nothing, # random placement of MSs with greedy user association
    "aux_assignment_params" => [
        "max_MSs_per_BS" => 1,

        "IA_infeasible_negative_inf_throughput" => false,

        "CoalitionFormationClustering_Individual:search_order" => :lexicographic,
        "CoalitionFormationClustering_Individual:stability_type" => :individual,
        "CoalitionFormationClustering_Individual:use_history" => true,
        "CoalitionFormationClustering_Individual:starting_point" => :grand,
    ],
    "aux_precoding_params" => [
        "initial_precoders" => "eigendirection",
        "stop_crit" => 1e-3,
        "max_iters" => 1000,
    ],
]
