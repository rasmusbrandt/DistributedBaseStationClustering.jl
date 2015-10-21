#!/usr/bin/env julia

##########################################################################
# timing-assignment.jl
#
# Timing for cluster assignment methods
##########################################################################

using DistributedBaseStationClustering, CoordinatedPrecoding
using Compat, JLD

##########################################################################
# General settings
srand(973472333)

##########################################################################
# Indoors network
simulation_params = @compat Dict(
    "I" => 8, "Kc" => 1, "N" => 2, "M" => 4, "d" => 1,
    "geography_size" => (1300.,1300.),
    "MS_serving_BS_distance" => Nullable{Float64}(),
    "Ntest" => 100,
    "assignment_methods" => [
        ExhaustiveSearchClustering,
        BranchAndBoundClustering,

        CoalitionFormationClustering_Group,
        CoalitionFormationClustering_Individual,

        GreedyClustering_Single,
        GreedyClustering_Multiple,

        Chen2014_LinearObj_ExhaustiveSearch,

        GrandCoalitionClustering,
        RandomClustering,
        NoClustering,
    ],
    "aux_network_params" => Dict(
        "num_coherence_symbols" => 2_700,
    ),
    "aux_assignment_params" => Dict(
        "clustering_type" => :spectrum_sharing,
        "apply_overhead_prelog" => false,
        "IA_infeasible_negative_inf_utility" => true,
        "replace_E1_utility_with_lower_bound" => false,

        "CoalitionFormationClustering_Group:max_merge_size" => 4,
        "CoalitionFormationClustering_Group:search_order" => :greedy,
        "CoalitionFormationClustering_Individual:search_budget" => 100,
        "CoalitionFormationClustering_Individual:search_order" => :greedy,
    ),
)
network =
    setup_random_large_scale_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        num_streams=simulation_params["d"],
        geography_size=simulation_params["geography_size"],
        MS_serving_BS_distance=simulation_params["MS_serving_BS_distance"])

timing(network, simulation_params, loop_over=:assignment_methods)
