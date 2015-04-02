#!/usr/bin/env julia

##########################################################################
# run_assignment-SNR.jl
#
# Performance as a function of transmit power, comparing different
# cluster assignment methods.
##########################################################################

include("src/IAClustering.jl")
using IAClustering, CoordinatedPrecoding
using HDF5, JLD

##########################################################################
# General settings
srand(8352339)
start_time = strftime("%Y%m%dT%H%M%S", time())

##########################################################################
# Indoors network
simulation_params = [
    "simulation_name" => "assignment-SNR_$(start_time)",
    "I" => 8, "Kc" => 1, "N" => 2, "M" => 2,
    "d" => 1,
    "Ndrops" => 10,
    "assignment_methods" => [
        # ExhaustiveSearchClustering,
        BranchAndBoundClustering,

        CoalitionFormationClustering_Group,
        CoalitionFormationClustering_Individual,

        Chen2014_ExhaustiveSearch,

        GrandCoalitionClustering,
        GreedyClustering,
        # RandomClustering,
        NoClustering,
    ],
    "aux_network_params" => [
        "no_coherence_symbols" => 1000,
    ],
    "aux_assignment_params" => [
        "clustering_type" => :spectrum_sharing,
        "apply_overhead_prelog" => false,
        "IA_infeasible_negative_inf_utility" => true,

        "CoalitionFormationClustering_Group:max_merge_size" => 4,
        "CoalitionFormationClustering_Group:search_order" => :greedy,
        "CoalitionFormationClustering_Individual:search_budget" => 100,
        "CoalitionFormationClustering_Individual:search_order" => :greedy,
    ],
    "independent_variable" => (set_transmit_powers_dBm!, -50:10:0),
]
network =
    setup_random_large_scale_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        no_streams=simulation_params["d"])
raw_results = simulate_assignment(network, simulation_params)

println("-- Saving $(simulation_params["simulation_name"]) results")
save("$(simulation_params["simulation_name"]).jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "raw_results", raw_results)
