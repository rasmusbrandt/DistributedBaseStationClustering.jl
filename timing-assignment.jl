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
    "I" => 8, "K" => 8, "N" => 2, "M" => 2,
    "d" => 1,
    "Ntest" => 100,
    "assignment_methods" => [
        ExhaustiveSearchClustering,
        BranchAndBoundClustering,

        CoalitionFormationClustering_Group,
        CoalitionFormationClustering_Individual,

        Chen2014_ExhaustiveSearch,

        GrandCoalitionClustering,
        GreedyClustering,
        RandomClustering,
        NoClustering,
    ],
    "aux_network_params" => [
        "no_coherence_symbols" => 1000,
    ],
    "aux_assignment_params" => [
        "clustering_type" => :spectrum_sharing,
        "apply_overhead_prelog" => false,
        "IA_infeasible_negative_inf_utility" => true,
        "replace_E1_utility_with_lower_bound" => false,

        "CoalitionFormationClustering_Group:max_merge_size" => 4,
        "CoalitionFormationClustering_Group:search_order" => :greedy,
        "CoalitionFormationClustering_Individual:search_budget" => 100,
        "CoalitionFormationClustering_Individual:search_order" => :greedy,
    ],
]
network =
    setup_random_large_scale_network(simulation_params["I"],
        simulation_params["K"], simulation_params["N"], simulation_params["M"],
        no_streams=simulation_params["d"])

timing(network, simulation_params, loop_over=:assignment_methods)
