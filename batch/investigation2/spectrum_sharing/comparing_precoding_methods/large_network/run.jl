#!/usr/bin/env julia

SRAND_SEED = 9825242

include(joinpath(dirname(@__FILE__), "../../../../../src/IAClustering.jl"))
using IAClustering, CoordinatedPrecoding
using HDF5, JLD

include(joinpath(dirname(@__FILE__), "../../../simulation_params.jl"))
include(joinpath(dirname(@__FILE__), "../../../simulation_params-large_network.jl"))
include(joinpath(dirname(@__FILE__), "../../../simulation_params-SNR.jl"))
include(joinpath(dirname(@__FILE__), "../../../simulation_params-precoding_methods.jl"))
include(joinpath(dirname(@__FILE__), "../../../plot_params-precoding_methods.jl"))

##########################################################################
# Plot setup
plot_params["axes"][:xlabel] = "MS speed [km/h]"

##########################################################################
# Simulation (both MinWLI and WMMSE)
start_time = strftime("%Y%m%dT%H%M%S", time())

srand(SRAND_SEED)
simulation_params["simulation_name"] = "precoding_methods-large_network-precoding_$(start_time)"
network =
    setup_random_large_scale_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        no_streams=simulation_params["d"],
        geography_size=simulation_params["geography_size"],
        MS_serving_BS_distance=simulation_params["MS_serving_BS_distance"])
raw_precoding_results, raw_assignment_results =
    simulate(network, simulation_params, loop_over=:precoding_methods)

println("-- Saving $(simulation_params["simulation_name"]) results")
save("$(simulation_params["simulation_name"]).jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "raw_precoding_results", raw_precoding_results,
     "raw_assignment_results", raw_assignment_results)

##########################################################################
# Generic plots
processed_results = postprocess(raw_precoding_results, simulation_params, plot_params)
plot(processed_results, simulation_params, plot_params)
