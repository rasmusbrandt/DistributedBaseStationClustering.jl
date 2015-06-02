v_kmh = 30 # km/h
v = v_kmh*(1e3/3600) # m/s

fd = v/(λ*Wc)
L = 1./(2*fd);

simulation_params["independent_variable"] = (set_average_SNRs_dB!, -10:10:50)
simulation_params["aux_network_params"] = [ "no_coherence_symbols" => L ]
