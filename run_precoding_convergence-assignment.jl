#!/usr/bin/env julia

##########################################################################
# run_precoding_convergence-assignment.jl
#
# Convergence, comparing different cluster assignment methods for the same
# precoding method.
##########################################################################

include("src/IAClustering.jl")
using IAClustering, CoordinatedPrecoding
using HDF5, JLD

##########################################################################
# Custom logging
Lumberjack.add_truck(Lumberjack.LumberjackTruck("debug.log", "debug"), "debug")

##########################################################################
# General settings
srand(83196723)
start_time = strftime("%Y%m%dT%H%M%S", time())

##########################################################################
# RandomLargeScaleNetwork
simulation_params = [
    "simulation_name" => "precoding_convergence-assignment_$(start_time)",
    "I" => 8, "Kc" => 1, "N" => 2, "M" => 2, "d" => 1,
    "Ndrops" => 10, "Nsim" => 20,
    "geography_length" => 250.,
    "MS_serving_BS_distance" => 50.,
    "assignment_methods" => [
        # ExhaustiveSearchClustering,
        BranchAndBoundClustering,

        CoalitionFormationClustering_Group,
        CoalitionFormationClustering_Individual,

        Chen2014_ExhaustiveSearch,
        # Peters2012_Heuristic,

        GrandCoalitionClustering,
        GreedyClustering,
        RandomClustering,
        NoClustering,
    ],
    "precoding_methods" => [
        RobustIntraclusterWMMSE,
    ],
    "aux_network_params" => [
        "no_coherence_symbols" => 1000,
    ],
    "aux_assignment_params" => [
        "clustering_type" => :spectrum_sharing,
        "apply_overhead_prelog" => true,
        "IA_infeasible_negative_inf_utility" => true,
        "replace_E1_utility_with_lower_bound" => false,

        "CoalitionFormationClustering_Group:max_merge_size" => 3,
        "CoalitionFormationClustering_Group:search_order" => :greedy,
        "CoalitionFormationClustering_Individual:search_budget" => 50,
        "CoalitionFormationClustering_Individual:search_order" => :greedy,
        "CoalitionFormationClustering_Individual:stability_type" => :contractual,
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
    setup_random_large_scale_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        no_streams=simulation_params["d"],
        geography_width=simulation_params["geography_length"],
        geography_height=simulation_params["geography_length"],
        MS_serving_BS_distance=simulation_params["MS_serving_BS_distance"])

raw_results =
    simulate_precoding_convergence(network, simulation_params, loop_over=:assignment_methods)

println("-- Saving $(simulation_params["simulation_name"]) results")
save("$(simulation_params["simulation_name"]).jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "raw_results", raw_results)
