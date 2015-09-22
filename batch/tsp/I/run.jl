#!/usr/bin/env julia

include(joinpath(dirname(@__FILE__), "../../../src/DistributedBSClustering.jl"))
using DistributedBSClustering, CoordinatedPrecoding
using Compat, JLD

# Parameters
include(joinpath(dirname(@__FILE__), "../simulation_params.jl"))
include(joinpath(dirname(@__FILE__), "../simulation_params-I.jl"))

# Ndrops are simulated below
simulation_params["Ndrops"] = 1
simulation_params["Nsim"] = 1

const Is = 1:32
const Ndrops = 1000

srand(725242)

results_searches = zeros(Float64, length(Is), Ndrops, 2)
results_throughputs = zeros(Float64, length(Is), Ndrops, 2)
for (idx, It) in enumerate(Is)
    network =
        setup_random_large_scale_network(It,
            simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
            num_streams=simulation_params["d"],
            geography_size=(sqrt(It*BS_density), sqrt(It*BS_density)),
            MS_serving_BS_distance=simulation_params["MS_serving_BS_distance"])

    for ii_Ndrop = 1:Ndrops
        _, raw_assignment_results =
            simulate(network, simulation_params, loop_over=:assignment_methods)
        results_searches[idx, ii_Ndrop, 1] = mean(raw_assignment_results[1]["CoalitionFormationClustering_AttachOrSupplant"]["num_searches"])
        results_searches[idx, ii_Ndrop, 2] = mean(raw_assignment_results[1]["CoalitionFormationClustering_Attach"]["num_searches"])
        results_throughputs[idx, ii_Ndrop, 1] = sum(raw_assignment_results[1]["CoalitionFormationClustering_AttachOrSupplant"]["throughputs"])
        results_throughputs[idx, ii_Ndrop, 2] = sum(raw_assignment_results[1]["CoalitionFormationClustering_Attach"]["throughputs"])
    end
end

println("-- Saving results")
save("I.jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "Is", Is,
     "results_searches", results_searches,
     "results_throughputs", results_throughputs)
