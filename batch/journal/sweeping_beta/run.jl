#!/usr/bin/env julia

include(joinpath(dirname(@__FILE__), "../../../src/IAClustering.jl"))
using CoordinatedPrecoding, IAClustering, JLD
Lumberjack.remove_truck("default")

include(joinpath(dirname(@__FILE__), "../simulation_params.jl"))
include(joinpath(dirname(@__FILE__), "../simulation_params-assignment_methods.jl"))
include(joinpath(dirname(@__FILE__), "../simulation_params-beta.jl"))
include(joinpath(dirname(@__FILE__), "../plot_params-assignment_methods.jl"))
include(joinpath(dirname(@__FILE__), "../plot_params-final.jl"))

##########################################################################
# Simulation
srand(927272)
simulation_params["simulation_name"] = "beta"
network =
    setup_random_large_scale_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        no_streams=simulation_params["d"],
        geography_size=simulation_params["geography_size"],
        MS_serving_BS_distance=simulation_params["MS_serving_BS_distance"])
raw_precoding_results, raw_assignment_results =
    simulate(network, simulation_params, loop_over=:assignment_methods)

##########################################################################
# Quick-and-dirty plots
processed_assignment_results = Dict{ASCIIString, Any}()
for p in (plot_params_longterm_sumrate, plot_params_longterm_avg_cluster_size, plot_params_longterm_num_sum_throughput_calculations, plot_params_longterm_num_searches)
    p["axes"][:xlabel] = "beta"
    tmp_processed_assignment_results = postprocess(raw_assignment_results, simulation_params, p)
    plot(tmp_processed_assignment_results, simulation_params, p)
end
plot_params_instantaneous_sumrate["axes"][:xlabel] = "beta"
tmp_processed_precoding_results = postprocess(raw_precoding_results, simulation_params, plot_params_instantaneous_sumrate)
plot(tmp_processed_precoding_results, simulation_params, plot_params_instantaneous_sumrate)

##########################################################################
# Save for publication plots
processed_assignment_results = postprocess(raw_assignment_results, simulation_params, plot_params_longterm_sumrate)
save("$(simulation_params["simulation_name"]).jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "processed_assignment_results", postprocess(raw_assignment_results, simulation_params, postprocess_params_assignment),
     "processed_precoding_results", postprocess(raw_precoding_results, simulation_params, postprocess_params_precoding))
