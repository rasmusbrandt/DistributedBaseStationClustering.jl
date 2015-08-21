#!/usr/bin/env julia

##########################################################################
# timing-assignment.jl
#
# Timing for cluster assignment methods
##########################################################################

include("src/IAClustering.jl")
using IAClustering, CoordinatedPrecoding
using Compat, JLD

##########################################################################
# General settings
srand(973472333)

##########################################################################
# Indoors network
simulation_params = [
    "I" => 8, "Kc" => 1, "N" => 2, "M" => 2,
    "d" => 1,
    "Ntest" => 100,
    "geography_size" => (1500.,1500.),
    "MS_serving_BS_distance" => Nullable(150.),
    "assignment_methods" => [
        ExhaustiveSearchClustering,
        BranchAndBoundClustering,

        CoalitionFormationClustering_AttachOrSupplant,
        CoalitionFormationClustering_Attach,

        GreedyClustering_Single,
        GreedyClustering_Multiple,

        Chen2014_LinearObj_ExhaustiveSearch,
        Chen2014_kmeans,
        Peters2012_Heuristic,

        GrandCoalitionClustering,
        RandomClustering,
        NoClustering,
    ],
    "aux_network_params" => [
        "num_coherence_symbols" => 2500,
        "beta_network_sdma" => 0.8,
    ],
    "aux_assignment_params" => [
        "BranchAndBoundClustering:max_abs_optimality_gap" => 0.,
        "BranchAndBoundClustering:E1_bound_in_rate_bound" => false,

        "CoalitionFormationClustering:search_order" => :random,
        "CoalitionFormationClustering:stability_type" => :individual,
        "CoalitionFormationClustering:search_budget" => 100,
        "CoalitionFormationClustering:use_history" => true,
        "CoalitionFormationClustering:starting_point" => :singletons,
    ],
]
network =
    setup_random_large_scale_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        num_streams=simulation_params["d"],
        geography_size=simulation_params["geography_size"],
        MS_serving_BS_distance=simulation_params["MS_serving_BS_distance"])

timing(network, simulation_params, loop_over=:assignment_methods)
