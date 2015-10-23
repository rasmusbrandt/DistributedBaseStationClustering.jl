#!/usr/bin/env julia

SRAND_SEED = 927272

include(joinpath(dirname(@__FILE__), "../../../../../src/DistributedBaseStationClustering.jl"))
using DistributedBaseStationClustering, CoordinatedPrecoding
using Compat, JLD

include(joinpath(dirname(@__FILE__), "../../../simulation_params.jl"))
include(joinpath(dirname(@__FILE__), "../../../simulation_params-assignment_methods.jl"))
include(joinpath(dirname(@__FILE__), "../../../simulation_params-small_network2.jl"))
include(joinpath(dirname(@__FILE__), "../../../simulation_params-SNR.jl"))
include(joinpath(dirname(@__FILE__), "../../../plot_params-assignment_methods.jl"))

##########################################################################
# Simulation
srand(SRAND_SEED)
network =
    setup_random_large_scale_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        num_streams=simulation_params["d"],
        geography_size=simulation_params["geography_size"],
        MS_serving_BS_distance=simulation_params["MS_serving_BS_distance"])
raw_precoding_results, raw_assignment_results =
    simulate(network, simulation_params, loop_over=:assignment_methods)

##########################################################################
# Generic plots
for p in (plot_params_longterm_sumrate, plot_params_longterm_num_utility_calculations, plot_params_longterm_num_clusters, plot_params_longterm_num_searches)
    processed_results = postprocess(raw_assignment_results, simulation_params, p)
    plot(processed_results, simulation_params, p)
end
for p in (plot_params_instantaneous_full_sumrate, plot_params_instantaneous_partial_sumrate, plot_params_instantaneous_LB_sumrate)
    processed_results = postprocess(raw_precoding_results, simulation_params, p)
    plot(processed_results, simulation_params, p)
end

##########################################################################
# Save for specialized plots
println("-- Saving $(simulation_params["simulation_name"]) results")
save("$(simulation_params["simulation_name"]).jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "raw_precoding_results", raw_precoding_results,
     "raw_assignment_results", raw_assignment_results)
