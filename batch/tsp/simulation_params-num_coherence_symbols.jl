const vs_kmh = linspace(0, 90, 46) # km/h
const vs = vs_kmh*(1e3/3600) # m/s

const fds = vs/(λ*Wc)
const Ls = 1./(2*fds)

simulation_params["independent_variable"] = ((n, v) -> set_aux_network_param!(n, v, "num_coherence_symbols"), Ls)
simulation_params["aux_independent_variables"] = [
    ((n, v) -> set_average_SNRs_dB!(n, v), [SNR_dB]),
    ((n, v) -> set_aux_network_param!(n, v, "beta_network_sdma"), [beta_network_sdma]),
]

push!(simulation_params["assignment_methods"], CoalitionFormationClustering_AttachOrSupplant_IgnoreIAFeasibility)
push!(simulation_params["assignment_methods"], CoalitionFormationClustering_Attach_IgnoreIAFeasibility)
