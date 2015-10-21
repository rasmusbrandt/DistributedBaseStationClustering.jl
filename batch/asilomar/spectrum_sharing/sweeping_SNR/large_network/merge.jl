#!/usr/bin/env julia

include(joinpath(dirname(@__FILE__), "../../../../../src/DistributedBaseStationClustering.jl"))
using DistributedBaseStationClustering, CoordinatedPrecoding
using Compat, JLD
using LaTeXStrings

include(joinpath(dirname(@__FILE__), "../../../simulation_params.jl"))
include(joinpath(dirname(@__FILE__), "../../../simulation_params-assignment_methods.jl"))
include(joinpath(dirname(@__FILE__), "../../../simulation_params-large_network.jl"))
include(joinpath(dirname(@__FILE__), "../../../simulation_params-SNR.jl"))
include(joinpath(dirname(@__FILE__), "../../../plot_params-assignment_methods.jl"))
include(joinpath(dirname(@__FILE__), "../../../plot_params-final.jl"))

##########################################################################
# Load data
sim_names = [
    "raw-large_network-run1.jld",
    "raw-large_network-run2.jld",
    "raw-large_network-run3.jld",
    "raw-large_network-run4.jld",
    "raw-large_network-run5.jld",
]

# Load first
println("Loading from $(sim_names[1])")
data = load(sim_names[1])
simulation_params = data["simulation_params"]
raw_assignment_results = data["raw_assignment_results"]
raw_precoding_results = data["raw_precoding_results"]

for sim_name in sim_names[2:end]
    println("Loading from $(sim_name)")
    data = load(sim_name)
    raw_assignment_results.simulation_results = cat(1, raw_assignment_results.simulation_results, data["raw_assignment_results"].simulation_results)
    raw_precoding_results.simulation_results = cat(1, raw_precoding_results.simulation_results, data["raw_precoding_results"].simulation_results)
end
data = []

##########################################################################
# Perform post processing
results_assignment, results_assignment_mean, results_assignment_var = postprocess(raw_assignment_results, simulation_params, postprocess_params_assignment)
results_precoding, results_precoding_mean, results_precoding_var = postprocess(raw_precoding_results, simulation_params, postprocess_params_precoding)

println("-- Saving merged results")
save("merged.jld",
     "simulation_params", simulation_params,
     "results_assignment", results_assignment,
     "results_assignment_mean", results_assignment_mean,
     "results_assignment_var", results_assignment_var,
     "results_precoding", results_precoding,
     "results_precoding_mean", results_precoding_mean,
     "results_precoding_var", results_precoding_var)
