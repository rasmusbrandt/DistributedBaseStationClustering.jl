fc = 2e9 # GHz
Wc = 300e3 # kHz

c = 300e6 # m/s
λ = c/fc # m

simulation_params = [
    "Ndrops" => 10, "Nsim" => 5,
    "geography_size" => (1500.,1500.),
    "MS_serving_BS_distance" => 150.,
    "aux_network_params" => [
        "beta_network_sdma" => 0.8,
    ],
    "aux_assignment_params" => Dict{ASCIIString, Any}(),
    "aux_precoding_params" => [
        "initial_precoders" => "eigendirection",
        "stop_crit" => 1e-2,
        "max_iters" => 1000,
    ],
]
