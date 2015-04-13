#!/usr/bin/env julia

include(joinpath(dirname(@__FILE__), "../../../../../src/IAClustering.jl"))
using IAClustering, CoordinatedPrecoding
using HDF5, JLD

include(joinpath(dirname(@__FILE__), "../../../simulation_params.jl"))
include(joinpath(dirname(@__FILE__), "../../../simulation_params-large_network.jl"))
include(joinpath(dirname(@__FILE__), "../../../simulation_params-no_coherence_symbols.jl"))
include(joinpath(dirname(@__FILE__), "../../../plot_params-assignment_methods.jl"))

##########################################################################
# Plot setup
for p in (plot_params_instantaneous_coord_sumrate, plot_params_instantaneous_noncoord_sumrate, plot_params_longterm_sumrate, plot_params_longterm_no_utility_calculations, plot_params_longterm_no_clusters)
    p["axes"][:xlabel] = "MS speed [km/h]"
    p["xvals"] = vs_kmh
end

##########################################################################
# Simulation (both MinWLI and WMMSE)
start_time = strftime("%Y%m%dT%H%M%S", time())

srand(725242)
simulation_params["simulation_name"] = "SNR-large_network-assignment-MinWLI_$(start_time)"
simulation_params["precoding_methods"] = [ RobustIntraclusterLeakageMinimization ]
network =
    setup_random_large_scale_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        no_streams=simulation_params["d"],
        geography_size=simulation_params["geography_size"],
        MS_serving_BS_distance=simulation_params["MS_serving_BS_distance"])
raw_precoding_results1, raw_assignment_results1 =
    simulate(network, simulation_params, loop_over=:assignment_methods)

processed_results = postprocess(raw_precoding_results1, simulation_params, plot_params_instantaneous_coord_sumrate)
plot(processed_results, simulation_params, plot_params_instantaneous_coord_sumrate)
processed_results = postprocess(raw_precoding_results1, simulation_params, plot_params_instantaneous_noncoord_sumrate)
plot(processed_results, simulation_params, plot_params_instantaneous_noncoord_sumrate)

println("-- Saving $(simulation_params["simulation_name"]) results")
save("$(simulation_params["simulation_name"]).jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "raw_precoding_results", raw_precoding_results1,
     "raw_assignment_results", raw_assignment_results1)

simulation_params["simulation_name"] = "SNR-large_network-assignment-WMMSE_$(start_time)"
simulation_params["precoding_methods"] = [ RobustIntraclusterWMMSE ]
network =
    setup_random_large_scale_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        no_streams=simulation_params["d"],
        geography_size=simulation_params["geography_size"],
        MS_serving_BS_distance=simulation_params["MS_serving_BS_distance"])
raw_precoding_results2, raw_assignment_results2 =
    simulate(network, simulation_params, loop_over=:assignment_methods)

processed_results = postprocess(raw_precoding_results2, simulation_params, plot_params_instantaneous_coord_sumrate)
plot(processed_results, simulation_params, plot_params_instantaneous_coord_sumrate)
processed_results = postprocess(raw_precoding_results2, simulation_params, plot_params_instantaneous_noncoord_sumrate)
plot(processed_results, simulation_params, plot_params_instantaneous_noncoord_sumrate)

println("-- Saving $(simulation_params["simulation_name"]) results")
save("$(simulation_params["simulation_name"]).jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "raw_precoding_results", raw_precoding_results2,
     "raw_assignment_results", raw_assignment_results2)

##########################################################################
# Generic plots
simulation_params["simulation_name"] = "SNR-large_network-assignment_$(start_time)"
for p in (plot_params_longterm_sumrate, plot_params_longterm_no_utility_calculations, plot_params_longterm_no_clusters)
    processed_results = postprocess(raw_assignment_results1, simulation_params, p)
    plot(processed_results, simulation_params, p)
end
