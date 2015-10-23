#!/usr/bin/env julia

using DistributedBaseStationClustering, CoordinatedPrecoding
using Compat, JLD

SEED = 7231936
srand(SEED)

simulation_params = @compat Dict(
    "simulation_name" => "movie-$SEED",
    "I" => 10, "Kc" => 1, "N" => 2, "M" => 3, "d" => 1,
    "Ndrops" => 1, "Nsim" => 1,
    "geography_size" => (1200.,1200.),
    "MS_serving_BS_distance" => Nullable{Float64}(),
    "assignment_methods" => [
        # BranchAndBoundClustering,
        CoalitionFormationClustering_Individual,

        GrandCoalitionClustering,
        NoClustering,
    ],
    "precoding_methods" => [ NoPrecoding, ],
    "aux_network_params" => Dict(
        "num_coherence_symbols" => 2_000,
    ),
    "aux_assignment_params" => Dict(
        "clustering_type" => :spectrum_sharing,
        "apply_overhead_prelog" => true,
        "IA_infeasible_negative_inf_utility" => false,
        "replace_E1_utility_with_lower_bound" => false,

        "BranchAndBoundClustering:bracket_E1" => false,

        "CoalitionFormationClustering_Individual:search_budget" => 100,
        "CoalitionFormationClustering_Individual:search_order" => :lexicographic,
        "CoalitionFormationClustering_Individual:stability_type" => :individual,
        "CoalitionFormationClustering_Individual:starting_point" => :grand,
        "CoalitionFormationClustering_Individual:use_history" => true,
        "CoalitionFormationClustering_Individual:prepare_movie" => true,
    ),
    "aux_precoding_params" => Dict("max_iters" => 1000),
    "independent_variable" => (set_transmit_powers_dBm!, [40.]),
)
network =
    setup_random_large_scale_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        num_streams=simulation_params["d"],
        geography_size=simulation_params["geography_size"],
        MS_serving_BS_distance=simulation_params["MS_serving_BS_distance"])

raw_precoding_results, raw_assignment_results =
    simulate(network, simulation_params, loop_over=:assignment_methods)

dir = simulation_params["simulation_name"]
output_filename = simulation_params["simulation_name"]
try
    mkdir(dir)
end
movie_state = raw_assignment_results[1]["CoalitionFormationClustering_Individual"]["movie_state"]
generate_movie(network, movie_state, dir, output_filename, 1.5)
