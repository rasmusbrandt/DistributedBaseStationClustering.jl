#!/usr/bin/env julia

include("../../../src/IAClustering.jl")
using IAClustering, CoordinatedPrecoding
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

        "BranchAndBoundClustering:max_abs_optimality_gap" => 0.,
        "BranchAndBoundClustering:max_rel_optimality_gap" => 0.,
        "BranchAndBoundClustering:E1_bound_in_rate_bound" => false,
        "BranchAndBoundClustering:store_fathomed_subtree_sizes" => true,
    ],
    "aux_precoding_params" => [
        "initial_precoders" => "eigendirection",
        "stop_crit" => 1e-2,
        "max_iters" => 1000,
    ],
    "independent_variable" => (set_average_SNRs_dB!, -10:5:40),
]

srand(725242)
network =
    setup_random_large_scale_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        num_streams=simulation_params["d"],
        geography_size=simulation_params["geography_size"],
        MS_serving_BS_distance=simulation_params["MS_serving_BS_distance"])

_, raw_assignment_results =
    simulate(network, simulation_params, loop_over=:assignment_methods)

println("-- Saving results")
save("SNR.jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "raw_assignment_results", raw_assignment_results)
