simulation_params["independent_variable"] = (set_average_SNRs_dB!, -10:50)
simulation_params["aux_independent_variables"] = [
    ((n, v) -> set_aux_network_param!(n, v, "num_coherence_symbols"), [num_coherence_symbols]),
    ((n, v) -> set_aux_network_param!(n, v, "beta_network_sdma"), [beta_network_sdma]),
]
