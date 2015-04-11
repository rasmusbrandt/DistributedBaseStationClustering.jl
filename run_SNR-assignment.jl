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
# General settings
srand(725242)
start_time = strftime("%Y%m%dT%H%M%S", time())

##########################################################################
# RandomLargeScaleNetwork
simulation_params = [
    "simulation_name" => "SNR-assignment_$(start_time)",
    "I" => 8, "K" => 8, "N" => 2, "M" => 2, "d" => 1,
    "Ndrops" => 10, "Nsim" => 20,
    "geography_size" => (250.,250.),
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
        "stop_crit" => 1e-3,
        "max_iters" => 1000,
    ],
    "independent_variable" => (set_transmit_powers_dBm!, -50:10:0),
]
network =
    setup_random_large_scale_network(simulation_params["I"],
        simulation_params["K"], simulation_params["N"], simulation_params["M"],
        no_streams=simulation_params["d"],
        geography_size=simulation_params["geography_size"])

raw_precoding_results, raw_assignment_results =
    simulate(network, simulation_params, loop_over=:assignment_methods)

println("-- Saving $(simulation_params["simulation_name"]) results")
save("$(simulation_params["simulation_name"]).jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "raw_precoding_results", raw_precoding_results,
     "raw_assignment_results", raw_assignment_results)
