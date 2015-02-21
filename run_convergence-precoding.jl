#!/usr/bin/env julia

##########################################################################
# run_convergence-precoding.jl
#
# Convergence, comparing different precoding methods for the same
# cluster assignment method.
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
srand(83196723)
start_time = strftime("%Y%m%dT%H%M%S", time())

##########################################################################
# Indoors network
simulation_params = [
    "simulation_name" => "convergence_$(start_time)-precoding",
    "I" => 6, "Kc" => 1, "N" => 2, "M" => 2,
    "d" => 1,
    "Ndrops" => 10, "Nsim" => 5,
    "assignment_methods" => [
        ExhaustiveSearchClustering,
    ],
    "precoding_methods" => [
        RobustIntraclusterWMMSE,
        NaiveIntraclusterWMMSE,
        RobustChen2014_MaxSINR,
        NaiveChen2014_MaxSINR,
        Shi2011_WMMSE,
        Eigenprecoding,
    ],
    "aux_precoding_params" => [
        "initial_precoders" => "eigendirection",
        "stop_crit" => 0.,
        "max_iters" => 20,
    ],
    "aux_independent_variables" => [
        (set_transmit_powers_dBm!, [-30, -10]),
    ]
]
network =
    setup_indoors_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        no_streams=simulation_params["d"])

raw_results = simulate_convergence(network, simulation_params, loop_over=:precoding_methods)

println("-- Saving $(simulation_params["simulation_name"]) results")
save("$(simulation_params["simulation_name"]).jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "raw_results", raw_results)
