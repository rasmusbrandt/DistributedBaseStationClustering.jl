#!/usr/bin/env julia

##########################################################################
# run_SNR-assignment.jl
#
# Performance as a function of transmit power, comparing different
# cluster assignment methods.
##########################################################################

include("src/IAClustering.jl")
using IAClustering, CoordinatedPrecoding
using HDF5, JLD

##########################################################################
# Set up logging
Lumberjack.configure(Lumberjack._lumber_mill.timber_trucks["console"]; mode = "warn")
Lumberjack.add_truck(Lumberjack.LumberjackTruck("debug.log", "debug"), "debug")

##########################################################################
# General settings
srand(973472333)
start_time = strftime("%Y%m%dT%H%M%S", time())

##########################################################################
# Indoors network
simulation_params = [
    "simulation_name" => "SNR_$(start_time)-assignment",
    "I" => 6, "Kc" => 1, "N" => 2, "M" => 2,
    "d" => 1,
    "Ndrops" => 10, "Nsim" => 5,
    "assignment_methods" => [
        GrandCoalitionClustering,
        ExhaustiveSearchClustering,
        Chen2014_LinearObjClustering_Exhaustive,
        # RandomClustering,
    ],
    "precoding_methods" => [
        RobustIntraclusterWMMSE,
    ],
    "aux_precoding_params" => [
        "initial_precoders" => "eigendirection",
        "stop_crit" => 1e-3,
        "max_iters" => 1000,
    ],
    "independent_variable" => (set_transmit_powers_dBm!, -80:10:0),
    # "aux_independent_variables" => [
    #     ((n, v) -> set_aux_precoding_param!(n, v, "max_iters"), [10, 50]),
    # ]
]
network =
    setup_indoors_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        no_streams=simulation_params["d"])
raw_results = simulate(network, simulation_params, loop_over=:assignment_methods)

println("-- Saving $(simulation_params["simulation_name"]) results")
save("$(simulation_params["simulation_name"]).jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "raw_results", raw_results)