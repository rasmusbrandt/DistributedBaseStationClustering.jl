#!/usr/bin/env julia

include("../../../src/IAClustering.jl")
using IAClustering, CoordinatedPrecoding
using Compat, JLD

# Parameters
include("../simulation_params.jl")
const Ndrops = 100
const Is = 1:16

simulation_params = [
    "Ndrops" => 1, "Nsim" => Nsim,
    "I" => 1, "Kc" => Kc,
    "M" => M, "N" => N, "d" => d,
    "geography_size" => (geography_width, geography_width),
    "MS_serving_BS_distance" => MS_serving_BS_distance,
    "assignment_methods" => [
        BranchAndBoundClustering,
        GreedyClustering_Multiple,
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
    "independent_variable" => (set_average_SNRs_dB!, [40]),
]

srand(725242)

results = zeros(Int, length(Is), Ndrops, 2)
for (idx, I) in enumerate(Is)
    network =
        setup_random_large_scale_network(I,
            simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
            num_streams=simulation_params["d"],
            geography_size=simulation_params["geography_size"],
            MS_serving_BS_distance=simulation_params["MS_serving_BS_distance"])

    for ii_Ndrop = 1:Ndrops
        _, raw_assignment_results =
            simulate(network, simulation_params, loop_over=:assignment_methods)
        results[idx, ii_Ndrop, 1] = raw_assignment_results[1]["BranchAndBoundClustering"]["num_sum_throughput_calculations"]
        results[idx, ii_Ndrop, 2] = raw_assignment_results[1]["GreedyClustering_Multiple"]["num_sum_throughput_calculations"]
    end
end

println("-- Saving results")
save("I.jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "Is", Is,
     "results", results)
