#!/usr/bin/env julia

include(joinpath(dirname(@__FILE__), "../../../src/IAClustering.jl"))
using IAClustering, CoordinatedPrecoding
using HDF5, JLD

include(joinpath(dirname(@__FILE__), "../simulation_params.jl"))
include(joinpath(dirname(@__FILE__), "../simulation_params-assignment_methods.jl"))
include(joinpath(dirname(@__FILE__), "../simulation_params-SNR.jl"))
include(joinpath(dirname(@__FILE__), "../plot_params-assignment_methods.jl"))

##########################################################################
# Simulation setup
SRAND_SEED = 927272; start_time = strftime("%Y%m%dT%H%M%S", time())

##########################################################################
# Nash
srand(SRAND_SEED)
simulation_params["simulation_name"] = "stability_concepts-nash_$(start_time)"
simulation_params["aux_assignment_params"]["CoalitionFormationClustering:stability_type"] = :nash
network =
    setup_random_large_scale_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        num_streams=simulation_params["d"],
        geography_size=simulation_params["geography_size"],
        MS_serving_BS_distance=simulation_params["MS_serving_BS_distance"])
raw_precoding_results, raw_assignment_results =
    simulate(network, simulation_params, loop_over=:assignment_methods)

for p in (plot_params_longterm_sumrate, plot_params_longterm_avg_cluster_size, plot_params_longterm_num_sum_throughput_calculations, plot_params_longterm_num_searches)
    p["axes"][:xlabel] = "SNR [dB]"
    processed_results = postprocess(raw_assignment_results, simulation_params, p)
    plot(processed_results, simulation_params, p)
end

plot_params_instantaneous_sumrate["axes"][:xlabel] = "SNR [dB]"
processed_results = postprocess(raw_precoding_results, simulation_params, plot_params_instantaneous_sumrate)
plot(processed_results, simulation_params, plot_params_instantaneous_sumrate)

##########################################################################
# Individual
srand(SRAND_SEED)
simulation_params["simulation_name"] = "stability_concepts-individual_$(start_time)"
simulation_params["aux_assignment_params"]["CoalitionFormationClustering:stability_type"] = :individual
network =
    setup_random_large_scale_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        num_streams=simulation_params["d"],
        geography_size=simulation_params["geography_size"],
        MS_serving_BS_distance=simulation_params["MS_serving_BS_distance"])
raw_precoding_results, raw_assignment_results =
    simulate(network, simulation_params, loop_over=:assignment_methods)

for p in (plot_params_longterm_sumrate, plot_params_longterm_avg_cluster_size, plot_params_longterm_num_sum_throughput_calculations, plot_params_longterm_num_searches)
    p["axes"][:xlabel] = "SNR [dB]"
    processed_results = postprocess(raw_assignment_results, simulation_params, p)
    plot(processed_results, simulation_params, p)
end

plot_params_instantaneous_sumrate["axes"][:xlabel] = "SNR [dB]"
processed_results = postprocess(raw_precoding_results, simulation_params, plot_params_instantaneous_sumrate)
plot(processed_results, simulation_params, plot_params_instantaneous_sumrate)

##########################################################################
# Swapee
srand(SRAND_SEED)
simulation_params["simulation_name"] = "stability_concepts-swapee_$(start_time)"
simulation_params["aux_assignment_params"]["CoalitionFormationClustering:stability_type"] = :swapee
network =
    setup_random_large_scale_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        num_streams=simulation_params["d"],
        geography_size=simulation_params["geography_size"],
        MS_serving_BS_distance=simulation_params["MS_serving_BS_distance"])
raw_precoding_results, raw_assignment_results =
    simulate(network, simulation_params, loop_over=:assignment_methods)

for p in (plot_params_longterm_sumrate, plot_params_longterm_avg_cluster_size, plot_params_longterm_num_sum_throughput_calculations, plot_params_longterm_num_searches)
    p["axes"][:xlabel] = "SNR [dB]"
    processed_results = postprocess(raw_assignment_results, simulation_params, p)
    plot(processed_results, simulation_params, p)
end

plot_params_instantaneous_sumrate["axes"][:xlabel] = "SNR [dB]"
processed_results = postprocess(raw_precoding_results, simulation_params, plot_params_instantaneous_sumrate)
plot(processed_results, simulation_params, plot_params_instantaneous_sumrate)

##########################################################################
# Contractual
srand(SRAND_SEED)
simulation_params["simulation_name"] = "stability_concepts-contractual_$(start_time)"
simulation_params["aux_assignment_params"]["CoalitionFormationClustering:stability_type"] = :contractual
network =
    setup_random_large_scale_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        num_streams=simulation_params["d"],
        geography_size=simulation_params["geography_size"],
        MS_serving_BS_distance=simulation_params["MS_serving_BS_distance"])
raw_precoding_results, raw_assignment_results =
    simulate(network, simulation_params, loop_over=:assignment_methods)

for p in (plot_params_longterm_sumrate, plot_params_longterm_avg_cluster_size, plot_params_longterm_num_sum_throughput_calculations, plot_params_longterm_num_searches)
    p["axes"][:xlabel] = "SNR [dB]"
    processed_results = postprocess(raw_assignment_results, simulation_params, p)
    plot(processed_results, simulation_params, p)
end

plot_params_instantaneous_sumrate["axes"][:xlabel] = "SNR [dB]"
processed_results = postprocess(raw_precoding_results, simulation_params, plot_params_instantaneous_sumrate)
plot(processed_results, simulation_params, plot_params_instantaneous_sumrate)
