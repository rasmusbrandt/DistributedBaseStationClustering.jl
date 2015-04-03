#!/usr/bin/env julia

##########################################################################
# run_SNR-precoding.jl
#
# Performance as a function of transmit power, comparing different
# precoding methods.
##########################################################################

include("src/IAClustering.jl")
using IAClustering, CoordinatedPrecoding
using HDF5, JLD

##########################################################################
# General settings
srand(973472333)
start_time = strftime("%Y%m%dT%H%M%S", time())

##########################################################################
# Indoors network
simulation_params = [
    "simulation_name" => "SNR-precoding_$(start_time)",
    "I" => 8, "Kc" => 1, "N" => 2, "M" => 2,
    "d" => 1,
    "Ndrops" => 10, "Nsim" => 20,
    "assignment_methods" => [
        BranchAndBoundClustering,
    ],
    "precoding_methods" => [
        RobustIntraclusterWMMSE,
        NaiveIntraclusterWMMSE,

        RobustIntraclusterLeakageMinimization,
        NaiveIntraclusterLeakageMinimization,

        RobustChen2014_MaxSINR,
        NaiveChen2014_MaxSINR,

        Shi2011_WMMSE,
        Eigenprecoding,
    ],
    "aux_network_params" => [
        "no_coherence_symbols" => 1000,
    ],
    "aux_assignment_params" => [
        "clustering_type" => :spectrum_sharing,
        "apply_overhead_prelog" => true,
        "IA_infeasible_negative_inf_utility" => true,
        "force_E1_utility_lower_bound" => true,
    ],
    "aux_precoding_params" => [
        "initial_precoders" => "eigendirection",
        "stop_crit" => 1e-3,
        "max_iters" => 1000,
    ],
    "independent_variable" => (set_transmit_powers_dBm!, -50:10:0),
]
network =
    setup_random_large_scale_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        no_streams=simulation_params["d"])
raw_precoding_results, raw_assignment_results =
    simulate(network, simulation_params, loop_over=:precoding_methods)

println("-- Saving $(simulation_params["simulation_name"]) results")
save("$(simulation_params["simulation_name"]).jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "raw_precoding_results", raw_precoding_results,
     "raw_assignment_results", raw_assignment_results)