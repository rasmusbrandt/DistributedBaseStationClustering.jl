v_kmh = 15 # km/h
v = v_kmh*(1e3/3600) # m/s

fd = v/(λ*Wc)
L = 1./(2*fd);

simulation_params["independent_variable"] = (set_transmit_powers_dBm!, -80:10:0)
simulation_params["aux_network_params"] = [ "num_coherence_symbols" => L ]
