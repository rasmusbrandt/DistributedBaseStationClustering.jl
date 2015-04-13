simulation_params["independent_variable"] = (set_average_SNRs_dB!, -10:10:40)

fc = 2.6e9 # GHz
Wc = 300e3 # kHz

v_kmh = 15 # km/h
v = v_kmh*(1e3/3600) # m/s

c = 300e6 # m/s
λ = c/fc # m

fd = v/(λ*Wc)
L = 1./(2*fd);

simulation_params["aux_network_params"] = [ "no_coherence_symbols" => L ]
