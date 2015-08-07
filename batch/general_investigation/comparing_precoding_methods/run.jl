#!/usr/bin/env julia

include(joinpath(dirname(@__FILE__), "../../../src/IAClustering.jl"))
using IAClustering, CoordinatedPrecoding
using HDF5, JLD

include(joinpath(dirname(@__FILE__), "../simulation_params.jl"))
include(joinpath(dirname(@__FILE__), "../simulation_params-precoding_methods.jl"))
include(joinpath(dirname(@__FILE__), "../simulation_params-SNR.jl"))
include(joinpath(dirname(@__FILE__), "../plot_params-precoding_methods.jl"))

##########################################################################
# Simulation
srand(927272); start_time = strftime("%Y%m%dT%H%M%S", time())
simulation_params["simulation_name"] = "precoding_methods_$(start_time)"
network =
    setup_random_large_scale_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        no_streams=simulation_params["d"],
        geography_size=simulation_params["geography_size"],
        MS_serving_BS_distance=simulation_params["MS_serving_BS_distance"])
raw_precoding_results, raw_assignment_results =
    simulate(network, simulation_params, loop_over=:precoding_methods)

##########################################################################
# Plots
plot_params["axes"][:xlabel] = "SNR [dB]"
processed_results = postprocess(raw_precoding_results, simulation_params, plot_params)
plot(processed_results, simulation_params, plot_params)
