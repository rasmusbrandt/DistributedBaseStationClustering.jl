simulation_params["independent_variable"] = (set_average_SNRs_dB!, -10:10:40)

fc = 2.6e9 # GHz
Wc = 270e3 # kHz

v_kmh = 15 # km/h
v = v_kmh*(1e3/3600) # m/s

c = 300e6 # m/s
λ = c/fc # m

fd = v/(λ*Wc)
L = 1./(2*fd);

simulation_params["aux_network_params"] = [ "no_coherence_symbols" => L ]

simulation_params["precoding_methods"] = [
    RobustIntraclusterWMMSE,
    NaiveIntraclusterWMMSE,

    RobustIntraclusterLeakageMinimization,
    NaiveIntraclusterLeakageMinimization,

    RobustChen2014_MaxSINR,
    NaiveChen2014_MaxSINR,

    Shi2011_WMMSE,
    Eigenprecoding,
]
