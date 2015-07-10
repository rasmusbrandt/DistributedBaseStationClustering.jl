# Default values
fc = 2e9 # GHz
Wc = 300e3 # kHz
c = 300e6 # m/s
λ = c/fc # m
v_kmh = 30 # km/h
v = v_kmh*(1e3/3600) # m/s
fd = v/(λ*Wc)
num_coherence_symbols = 1/(2*fd)

beta_network_sdma = 0.8

SNR_dB = 20

simulation_params = [
    "Ndrops" => 10, "Nsim" => 5,
    "geography_size" => (1500.,1500.),
    "MS_serving_BS_distance" => 150.,
    "aux_network_params" => Dict{ASCIIString, Any}(),
    "aux_assignment_params" => Dict{ASCIIString, Any}(),
    "aux_precoding_params" => [
        "initial_precoders" => "eigendirection",
        "stop_crit" => 1e-2,
        "max_iters" => 1000,
    ],
]
