fc = 2.6e9 # GHz
Wc = 300e3 # kHz

vs_kmh = linspace(1, 100, 10) # km/h
vs = vs_kmh*(1e3/3600) # m/s

c = 300e6 # m/s
λ = c/fc # m

fds = vs/(λ*Wc)
Ls = 1./(2*fds);

simulation_params["independent_variable"] = ((n, v) -> set_aux_network_param!(n, v, "no_coherence_symbols"), Ls)
simulation_params["aux_independent_variables"] = [ ((n, v) -> set_average_SNRs_dB!(n, v), [20]) ]
