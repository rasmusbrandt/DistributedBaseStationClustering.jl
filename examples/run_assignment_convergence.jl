#!/usr/bin/env julia

##########################################################################
# run_assignment_convergence-assignment.jl
#
# Convergence of the assignment methods.
##########################################################################

using CoordinatedPrecoding, DistributedBaseStationClustering
using Compat, JLD

##########################################################################
# Custom logging
Lumberjack.add_truck(Lumberjack.LumberjackTruck("debug.log", "debug"), "debug")

##########################################################################
# General settings
srand(83196723)
start_time = strftime("%Y%m%dT%H%M%S", time())

##########################################################################
# RandomLargeScaleNetwork
simulation_params = @compat Dict(
    "simulation_name" => "assignment_convergence_$(start_time)",
    "I" => 12, "Kc" => 2, "N" => 2, "M" => 8, "d" => 1,
    "Ndrops" => 1,
    "geography_size" => (1500.,1500.),
    "MS_serving_BS_distance" => Nullable(150.),
    "assignment_methods" => [ BranchAndBoundClustering, ],
    "aux_network_params" => Dict(
        "num_coherence_symbols" => 2700,
        "beta_network_sdma" => 0.5,
    ),
    "aux_assignment_params" => Dict(
        "BranchAndBoundClustering:branching_rule" => :bfs,
        "BranchAndBoundClustering:max_abs_optimality_gap" => 0.,
        "BranchAndBoundClustering:max_rel_optimality_gap" => 0.,
        "BranchAndBoundClustering:E1_bound_in_rate_bound" => false,
        "BranchAndBoundClustering:store_evolution" => true,
    ),
)
network =
    setup_random_large_scale_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        num_streams=simulation_params["d"],
        geography_size=simulation_params["geography_size"],
        MS_serving_BS_distance=simulation_params["MS_serving_BS_distance"])

raw_results =
    simulate_assignment_convergence(network, simulation_params)

println("-- Saving $(simulation_params["simulation_name"]) results")
save("$(simulation_params["simulation_name"]).jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "raw_results", raw_results)
