vs_kmh = linspace(1, 100, 9) # km/h
vs = vs_kmh*(1e3/3600) # m/s

fds = vs/(λ*Wc)
Ls = 1./(2*fds);

transmit_power_dBm_sim = 25 - 10*log10(600)

simulation_params["independent_variable"] = ((n, v) -> set_aux_network_param!(n, v, "num_coherence_symbols"), Ls)
simulation_params["aux_independent_variables"] = [ ((n, v) -> set_transmit_powers_dBm!(n, v), [transmit_power_dBm_sim]) ]

simulation_params["aux_assignment_params"]["CoalitionFormationClustering_Individual:search_budget"] = 10
