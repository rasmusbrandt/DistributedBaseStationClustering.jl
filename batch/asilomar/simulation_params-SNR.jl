v_kmh = 30 # km/h
v = v_kmh*(1e3/3600) # m/s

fd = v/(λ*Wc)
L = 1./(2*fd);

transmit_powers_dBm = -10:5:60 # macro BS typicall has output power 46 dBm, pico has 20-30 dB less
transmit_powers_dBm_sim = transmit_powers_dBm - 10*log10(600) # simulating one 15 kHz subcarrier, of which there are 600 in a 10 MHz system

simulation_params["independent_variable"] = (set_transmit_powers_dBm!, transmit_powers_dBm_sim)
simulation_params["aux_independent_variables"] = [ ((n, v) -> set_aux_assignment_param!(n, v, "CoalitionFormationClustering_Individual:search_budget"), [2, 10]) ]
simulation_params["aux_network_params"] = [ "num_coherence_symbols" => L ]
