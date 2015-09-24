#!/usr/bin/env julia

using CoordinatedPrecoding, DistributedBaseStationClustering
using Compat, JLD

# Parameters
include("../simulation_params.jl")
const Ndrops = 10

simulation_params = [
    "Ndrops" => Ndrops, "Nsim" => Nsim,
    "I" => I, "Kc" => Kc,
    "M" => M, "N" => N, "d" => d,
    "geography_size" => (geography_width, geography_width),
    "MS_serving_BS_distance" => MS_serving_BS_distance,
    "assignment_methods" => [
        # ExhaustiveSearchClustering,
        BranchAndBoundClustering,

        CoalitionFormationClustering_AttachOrSupplant,

        GreedyClustering_Single,
        GreedyClustering_Multiple,

        GrandCoalitionClustering,
        NoClustering,
    ],
    "precoding_methods" => [ NoPrecoding ],
    "aux_network_params" => [
        "num_coherence_symbols" => num_coherence_symbols,
        "beta_network_sdma" => beta_network_sdma,
    ],
    "aux_assignment_params" => [
        "max_num_MSs_per_BS" => Kc,

        "BranchAndBoundClustering:branching_rule" => :dfs,
        "BranchAndBoundClustering:max_abs_optimality_gap" => 0.,
        "BranchAndBoundClustering:max_rel_optimality_gap" => 0.,
        "BranchAndBoundClustering:E1_bound_in_rate_bound" => false,
        "BranchAndBoundClustering:store_evolution" => false,
    ],
    "aux_precoding_params" => [
        "initial_precoders" => "eigendirection",
        "stop_crit" => 1e-2,
        "max_iters" => 1000,
    ],
    "independent_variable" => (set_average_SNRs_dB!, -10:5:40),
]
