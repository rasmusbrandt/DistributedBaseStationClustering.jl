simulation_params["independent_variable"] = ((n, v) -> set_aux_network_param!(n, v, "beta_network_sdma"), [1e-5, 0.2:0.2:1])
simulation_params["aux_independent_variables"] = [
    ((n, v) -> set_average_SNRs_dB!(n, v), [SNR_dB]),
    ((n, v) -> set_aux_network_param!(n, v, "num_coherence_symbols"), [num_coherence_symbols]),
]
