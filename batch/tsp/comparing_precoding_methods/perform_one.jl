#!/usr/bin/env julia

using CoordinatedPrecoding, DistributedBaseStationClustering
using Compat, JLD, ArgParse
Lumberjack.remove_truck("default")

include(joinpath(dirname(@__FILE__), "../simulation_params.jl"))
include(joinpath(dirname(@__FILE__), "../simulation_params-precoding_methods.jl"))
include(joinpath(dirname(@__FILE__), "../simulation_params-SNR.jl"))
include(joinpath(dirname(@__FILE__), "../plot_params-precoding_methods.jl"))

function perform_one()
    # Get seed
    s = ArgParseSettings()
    @add_arg_table s begin
        "--seed", "-s"
            help = "RNG seed"
            arg_type = Int
            required = true
    end
    args = parse_args(s)
    seed = args["seed"]

    # Simulation
    srand(seed)
    simulation_params["simulation_name"] = "precoding_methods-seed$seed"
    network =
        setup_random_large_scale_network(simulation_params["I"],
            simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
            num_streams=simulation_params["d"],
            geography_size=simulation_params["geography_size"],
            MS_serving_BS_distance=simulation_params["MS_serving_BS_distance"])
    raw_precoding_results, raw_assignment_results =
        simulate(network, simulation_params, loop_over=:precoding_methods)

    # Quick-and-dirty plots
    plot_params["axes"][:xlabel] = "SNR [dB]"
    tmp_processed_precoding_results = postprocess(raw_precoding_results, simulation_params, plot_params)
    plot(tmp_processed_precoding_results, simulation_params, plot_params)

    # Save for publication plots
    save("$(simulation_params["simulation_name"]).jld",
         "simulation_params", clean_simulation_params_for_jld(simulation_params),
         "raw_assignment_results", raw_assignment_results,
         "raw_precoding_results", raw_precoding_results)
end

perform_one()
