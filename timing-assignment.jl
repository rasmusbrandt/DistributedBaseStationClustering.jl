#!/usr/bin/env julia

##########################################################################
# timing-assignment.jl
#
# Timing for cluster assignment methods
##########################################################################

include("src/IAClustering.jl")
using IAClustering, CoordinatedPrecoding
using HDF5, JLD

##########################################################################
# General settings
srand(973472333)

##########################################################################
# Indoors network
simulation_params = [
    "I" => 12, "Kc" => 1, "N" => 2, "M" => 2,
    "d" => 1,
    "Ntest" => 10,
    "assignment_methods" => [
        Chen2014_LinearObjClustering_ExhaustiveSearch,
        ExhaustiveSearchClustering,
        BranchAndBoundClustering,

        GrandCoalitionClustering,
        RandomClustering,
        NoClustering,
    ],
    "precoding_methods" => [
        RobustIntraclusterWMMSE # pseudo
    ],
    "aux_precoding_params" => [
        "initial_precoders" => "eigendirection",
        "stop_crit" => 0.,
        "max_iters" => 100,
    ],
]
network =
    setup_indoors_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        no_streams=simulation_params["d"])

timing(network, simulation_params, loop_over=:assignment_methods)
