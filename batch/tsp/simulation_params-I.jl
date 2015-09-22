simulation_params["assignment_methods"] = [
    CoalitionFormationClustering_AttachOrSupplant,
    CoalitionFormationClustering_Attach,
]
simulation_params["precoding_methods"] = [ NoPrecoding, ]

simulation_params["independent_variable"] = (set_average_SNRs_dB!, [SNR_dB])
simulation_params["aux_independent_variables"] = [
    ((n, v) -> set_aux_network_param!(n, v, "num_coherence_symbols"), [10*num_coherence_symbols]),
    ((n, v) -> set_aux_network_param!(n, v, "beta_network_sdma"), [beta_network_sdma]),
]