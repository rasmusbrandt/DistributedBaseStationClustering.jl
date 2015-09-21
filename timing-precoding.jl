#!/usr/bin/env julia

##########################################################################
# timing-precoding.jl
#
# Timing for cluster precoding methods
##########################################################################

include("src/LongtermIAClustering.jl")
using LongtermIAClustering, CoordinatedPrecoding
using Compat, JLD

##########################################################################
# General settings
srand(973472333)

##########################################################################
# Indoors network
simulation_params = @Compat.Dict(
    "I" => 12, "Kc" => 2, "N" => 2, "M" => 8, "d" => 1,
    "Ndrops" => 10, "Nsim" => 5,
    "geography_size" => (1500.,1500.),
    "MS_serving_BS_distance" => Nullable(150.),
    "Ntest" => 100,
    "assignment_methods" => [ BranchAndBoundClustering ],
    "precoding_methods" => [
        RobustIntraclusterWMMSE,
        NaiveIntraclusterWMMSE,

        RobustIntraclusterLeakageMinimization,
        NaiveIntraclusterLeakageMinimization,

        RobustChen2014_MaxSINR,
        NaiveChen2014_MaxSINR,
    ],
    "aux_network_params" => @Compat.Dict(
        "num_coherence_symbols" => 2500,
        "beta_network_sdma" => 0.8,
    ),
    "aux_assignment_params" => @Compat.Dict(
        "BranchAndBoundClustering:branching_rule" => :dfs,
        "BranchAndBoundClustering:max_abs_optimality_gap" => 0.,
        "BranchAndBoundClustering:max_rel_optimality_gap" => 0.,
        "BranchAndBoundClustering:E1_bound_in_rate_bound" => false,
        "BranchAndBoundClustering:store_evolution" => false,
    ),
    "aux_precoding_params" => @Compat.Dict(
        "initial_precoders" => "eigendirection",
        "stop_crit" => 0.,
        "max_iters" => 100,
    ),
)
network =
    setup_random_large_scale_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        num_streams=simulation_params["d"],
        geography_size=simulation_params["geography_size"],
        MS_serving_BS_distance=simulation_params["MS_serving_BS_distance"])

timing(network, simulation_params, loop_over=:precoding_methods)
