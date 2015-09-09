#!/usr/bin/env julia

##########################################################################
# run_SNR-precoding.jl
#
# Performance as a function of transmit power, comparing different
# precoding methods.
##########################################################################

include("src/IAClustering.jl")
using IAClustering, CoordinatedPrecoding
using Compat, JLD

##########################################################################
# General settings
srand(973472333)
start_time = strftime("%Y%m%dT%H%M%S", time())

##########################################################################
# RandomLargeScaleNetwork
simulation_params = @Compat.Dict(
    "simulation_name" => "SNR-precoding_$(start_time)",
    "I" => 12, "Kc" => 2, "N" => 2, "M" => 8, "d" => 1,
    "Ndrops" => 10, "Nsim" => 5,
    "geography_size" => (1500.,1500.),
    "MS_serving_BS_distance" => Nullable(150.),
    "assignment_methods" => [ BranchAndBoundClustering, ],
    "precoding_methods" => [
        RobustIntraclusterWMMSE,
        NaiveIntraclusterWMMSE,

        RobustChen2014_MaxSINR,
        NaiveChen2014_MaxSINR,
    ],
    "aux_network_params" => @Compat.Dict(
        "num_coherence_symbols" => 2500,
        "beta_network_sdma" => 0.8,
    ),
    "aux_assignment_params" => @Compat.Dict(
        "BranchAndBoundClustering:max_abs_optimality_gap" => 0.,
        "BranchAndBoundClustering:max_rel_optimality_gap" => 0.,
        "BranchAndBoundClustering:E1_bound_in_rate_bound" => false,
    ),
    "aux_precoding_params" => @Compat.Dict(
        "initial_precoders" => "eigendirection",
        "stop_crit" => 1e-2,
        "max_iters" => 1000,
    ),
    "independent_variable" => (set_average_SNRs_dB!, -10:5:50),
)
network =
    setup_random_large_scale_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        num_streams=simulation_params["d"],
        geography_size=simulation_params["geography_size"],
        MS_serving_BS_distance=simulation_params["MS_serving_BS_distance"])

raw_precoding_results, raw_assignment_results =
    simulate(network, simulation_params, loop_over=:precoding_methods)

println("-- Saving $(simulation_params["simulation_name"]) results")
save("$(simulation_params["simulation_name"]).jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "raw_precoding_results", raw_precoding_results,
     "raw_assignment_results", raw_assignment_results)
