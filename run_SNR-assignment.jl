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
    "I" => 8, "Kc" => 1, "N" => 2, "M" => 2, "d" => 1,
    "Ndrops" => 10, "Nsim" => 5,
    "geography_size" => (1500.,1500.),
    "MS_serving_BS_distance" => 150.,
    "assignment_methods" => [
        # ExhaustiveSearchClustering,
        BranchAndBoundClustering,

        CoalitionFormationClustering_Swap,
        CoalitionFormationClustering_Individual,
        # CoalitionFormationClustering_Group,

        GreedyClustering_Single,
        GreedyClustering_Multiple,

        # Chen2014_LinearObj_ExhaustiveSearch,
        Chen2014_kmeans,
        # Peters2012_Heuristic,

        GrandCoalitionClustering,
        RandomClustering,
        NoClustering,
    ],
    "precoding_methods" => [ RobustIntraclusterWMMSE, ],
    "aux_network_params" => [
        "num_coherence_symbols" => 2500,
        "beta_network_sdma" => 0.8,
    ],
    "aux_assignment_params" => [
        "BranchAndBoundClustering:max_abs_optimality_gap" => 0.,
        "BranchAndBoundClustering:E1_bound_in_rate_bound" => false,

        "CoalitionFormationClustering_Swap:search_order" => :random,
        "CoalitionFormationClustering_Swap:stability_type" => :individual,
        "CoalitionFormationClustering_Swap:search_budget" => 100,
        "CoalitionFormationClustering_Swap:use_history" => true,
        "CoalitionFormationClustering_Swap:starting_point" => :singletons,

        "CoalitionFormationClustering_Individual:search_order" => :random,
        "CoalitionFormationClustering_Individual:stability_type" => :individual,
        "CoalitionFormationClustering_Individual:search_budget" => 100,
        "CoalitionFormationClustering_Individual:use_history" => true,
        "CoalitionFormationClustering_Individual:starting_point" => :singletons,

        "CoalitionFormationClustering_Group:max_num_merging_coalitions" => 3,
        "CoalitionFormationClustering_Group:search_order" => :greedy,
    ],
    "aux_precoding_params" => [
        "initial_precoders" => "eigendirection",
        "stop_crit" => 1e-2,
        "max_iters" => 1000,
    ],
    "independent_variable" => (set_average_SNRs_dB!, -10:5:50),
]
network =
    setup_random_large_scale_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        no_streams=simulation_params["d"],
        geography_size=simulation_params["geography_size"],
        MS_serving_BS_distance=simulation_params["MS_serving_BS_distance"])

raw_precoding_results, raw_assignment_results =
    simulate(network, simulation_params, loop_over=:assignment_methods)

println("-- Saving $(simulation_params["simulation_name"]) results")
save("$(simulation_params["simulation_name"]).jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "raw_precoding_results", raw_precoding_results,
     "raw_assignment_results", raw_assignment_results)
