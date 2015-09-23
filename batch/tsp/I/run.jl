#!/usr/bin/env julia

include(joinpath(dirname(@__FILE__), "../../../src/DistributedBSClustering.jl"))
using DistributedBSClustering, CoordinatedPrecoding
using Compat, JLD

# Parameters
include(joinpath(dirname(@__FILE__), "../simulation_params.jl"))
include(joinpath(dirname(@__FILE__), "../simulation_params-I.jl"))

srand(7315242)

const Is = 2:4:50
simulation_params["Ndrops"] *= 10 # no run.sh

results_searches = zeros(Float64, length(Is), simulation_params["Ndrops"], 2)
results_throughputs = zeros(Float64, length(Is), simulation_params["Ndrops"], 2)
for (idx, It) in enumerate(Is)
    network =
        setup_random_large_scale_network(It,
            simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
            num_streams=simulation_params["d"],
            geography_size=(sqrt(It*BS_density), sqrt(It*BS_density)),
            MS_serving_BS_distance=simulation_params["MS_serving_BS_distance"])

    _, raw_assignment_results =
        simulate(network, simulation_params, loop_over=:assignment_methods)
    results_searches[idx, :, 1] = [ mean(raw_assignment_results[ii_Ndrop]["CoalitionFormationClustering_AttachOrSupplant"]["num_searches"]) for ii_Ndrop = 1:simulation_params["Ndrops"] ]
    results_searches[idx, :, 2] = [ mean(raw_assignment_results[ii_Ndrop]["CoalitionFormationClustering_Attach"]["num_searches"]) for ii_Ndrop = 1:simulation_params["Ndrops"] ]
    results_throughputs[idx, :, 1] = [ sum(raw_assignment_results[ii_Ndrop]["CoalitionFormationClustering_AttachOrSupplant"]["throughputs"]) for ii_Ndrop = 1:simulation_params["Ndrops"] ]
    results_throughputs[idx, :, 2] = [ sum(raw_assignment_results[ii_Ndrop]["CoalitionFormationClustering_Attach"]["throughputs"]) for ii_Ndrop = 1:simulation_params["Ndrops"] ]
end

println("-- Saving results")
save("I.jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "Is", Is,
     "results_searches", results_searches,
     "results_throughputs", results_throughputs)
