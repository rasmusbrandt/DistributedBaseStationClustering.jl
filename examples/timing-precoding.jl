#!/usr/bin/env julia

##########################################################################
# timing-precoding.jl
#
# Timing for cluster precoding methods
##########################################################################

using DistributedBaseStationClustering, CoordinatedPrecoding
using Compat, JLD

##########################################################################
# General settings
srand(973472333)

##########################################################################
# Indoors network
simulation_params = @compat Dict(
    "I" => 8, "Kc" => 1, "N" => 2, "M" => 2, "d" => 1,
    "geography_size" => (1300.,1300.),
    "MS_serving_BS_distance" => Nullable{Float64}(),
    "Ntest" => 100,
    "assignment_methods" => [ GrandCoalitionClustering ],
    "precoding_methods" => [
        RobustIntraclusterWMMSE,
        NaiveIntraclusterWMMSE,

        RobustIntraclusterLeakageMinimization,
        NaiveIntraclusterLeakageMinimization,

        RobustChen2014_MaxSINR,
        NaiveChen2014_MaxSINR,

        Shi2011_WMMSE,
        Eigenprecoding,
    ],
    "aux_network_params" => Dict(
        "num_coherence_symbols" => 2_700,
    ),
    "aux_assignment_params" => Dict(
        "clustering_type" => :spectrum_sharing,
        "apply_overhead_prelog" => false,
        "IA_infeasible_negative_inf_utility" => true,
        "replace_E1_utility_with_lower_bound" => false,
    ),
    "aux_precoding_params" => Dict(
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
