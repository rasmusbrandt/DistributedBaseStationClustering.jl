#!/usr/bin/env julia

using CoordinatedPrecoding, DistributedBaseStationClustering
using Compat, JLD

include(joinpath(dirname(@__FILE__), "../plot_params-precoding_methods.jl"))

##########################################################################
# Load data
sim_name = "precoding_methods"
seeds = split(readchomp("../seeds.txt"), '\n'); Nseeds = length(seeds)
file_names = [ "$(sim_name)-seed$n.jld" for n in seeds ]

# Load first
println("Initial loading from $(file_names[1])")
data = load(file_names[1])
simulation_params = data["simulation_params"]; Ndrops = simulation_params["Ndrops"]
assn_siz = size(data["raw_assignment_results"].simulation_results)
prec_siz = size(data["raw_precoding_results"].simulation_results)

raw_assignment_results = CoordinatedPrecoding.MultipleSimulationResults(Nseeds*assn_siz[1], assn_siz[2:end]...)
raw_precoding_results = CoordinatedPrecoding.MultipleSimulationResults(Nseeds*assn_siz[1], prec_siz[2:end]...)

for (file_idx, file_name) in enumerate(file_names)
    println("Loading from $(file_name)")
    data = load(file_name)
    raw_precoding_results.simulation_results[(file_idx - 1)*Ndrops + 1:file_idx*Ndrops, :, :, :] = data["raw_precoding_results"].simulation_results
    data = []
end

##########################################################################
# Save merged and post processed data
simulation_params["Ndrops"] *= Nseeds
println("-- Saving merged results")
save("$(sim_name).jld",
     "simulation_params", simulation_params,
     "processed_precoding_results", postprocess(raw_precoding_results, simulation_params, plot_params))
