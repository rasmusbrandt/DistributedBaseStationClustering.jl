#!/usr/bin/env julia

include(joinpath(dirname(@__FILE__), "../../../src/IAClustering.jl"))
using CoordinatedPrecoding, IAClustering
using Compat, JLD, ArgParse

include(joinpath(dirname(@__FILE__), "../plot_params-final.jl"))

##########################################################################
# File names
s = ArgParseSettings()
@add_arg_table s begin
    "file_names"
        help = "file names with results"
        required = true
        nargs = '+'
end
parsed_args = parse_args(s)

##########################################################################
# Load data
sim_name = "precoding_methods"
file_names = parsed_args["file_names"]

# Load first
println("Loading from $(file_names[1])")
data = load(file_names[1])
simulation_params = data["simulation_params"]
raw_precoding_results = data["raw_precoding_results"]

for file_name in file_names[2:end]
    println("Loading from $(file_name)")
    data = load(file_name)
    simulation_params["Ndrops"] += data["simulation_params"]["Ndrops"]
    raw_precoding_results.simulation_results = cat(1, raw_precoding_results.simulation_results, data["raw_precoding_results"].simulation_results)
end
data = []

##########################################################################
# Save merged and post processed data
println("-- Saving merged results")
save("$(sim_name).jld",
     "simulation_params", simulation_params,
     "processed_precoding_results", postprocess(raw_precoding_results, simulation_params, postprocess_params_precoding2))
