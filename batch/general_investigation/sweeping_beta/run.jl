#!/usr/bin/env julia

include(joinpath(dirname(@__FILE__), "../../../src/IAClustering.jl"))
using IAClustering, CoordinatedPrecoding
using HDF5, JLD

include(joinpath(dirname(@__FILE__), "../simulation_params.jl"))
include(joinpath(dirname(@__FILE__), "../simulation_params-assignment_methods.jl"))
include(joinpath(dirname(@__FILE__), "../simulation_params-beta.jl"))
include(joinpath(dirname(@__FILE__), "../plot_params-assignment_methods.jl"))

##########################################################################
# Simulation
srand(927272); start_time = strftime("%Y%m%dT%H%M%S", time())
simulation_params["simulation_name"] = "beta_$(start_time)"
network =
    setup_random_large_scale_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        no_streams=simulation_params["d"],
        geography_size=simulation_params["geography_size"],
        MS_serving_BS_distance=simulation_params["MS_serving_BS_distance"])
raw_precoding_results, raw_assignment_results =
    simulate(network, simulation_params, loop_over=:assignment_methods)

##########################################################################
# Plot
for p in (plot_params_longterm_sumrate, plot_params_longterm_avg_cluster_size, plot_params_longterm_num_sum_throughput_calculations)
    p["axes"][:xlabel] = "beta"
    processed_results = postprocess(raw_assignment_results, simulation_params, p)
    plot(processed_results, simulation_params, p)
end

plot_params_instantaneous_sumrate["axes"][:xlabel] = "beta"
processed_results = postprocess(raw_precoding_results, simulation_params, plot_params_instantaneous_sumrate)
plot(processed_results, simulation_params, plot_params_instantaneous_sumrate)
