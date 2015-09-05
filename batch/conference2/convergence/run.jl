#!/usr/bin/env julia

include("../../../src/IAClustering.jl")
using IAClustering, CoordinatedPrecoding
using Compat, JLD

# Parameters
const fc = 2e9 # GHz
const Wc = 300e3 # kHz
const c = 300e6 # m/s
const λ = c/fc # m
const v_kmh = 30 # km/h
const v = v_kmh*(1e3/3600) # m/s
const fd = v/(λ*Wc)
const num_coherence_symbols = 1/(2*fd)
const beta_network_sdma = 0.5

const design_ISD = 500
const BS_density = sqrt(3)/2*design_ISD^2 # BSs per m^2, for hexagonal cells
const I = 16; const Kc = 2
const M = 8; const N = 2; const d = 1
const geography_width = sqrt(I*BS_density)
const MS_serving_BS_distance = Nullable(150.) # geography_width/10. = 161.1854897735313

const SNR_dB = 20

const Ndrops = 1
const Nsim = 1

simulation_params = [
    "Ndrops" => Ndrops, "Nsim" => Nsim,
    "I" => I, "Kc" => Kc,
    "M" => M, "N" => N, "d" => d,
    "geography_size" => (geography_width, geography_width),
    "MS_serving_BS_distance" => MS_serving_BS_distance,
    "assignment_methods" => [ BranchAndBoundClustering, ],
    "aux_network_params" => [
        "num_coherence_symbols" => num_coherence_symbols,
        "beta_network_sdma" => beta_network_sdma,
    ],
    "aux_assignment_params" => [
        "max_num_MSs_per_BS" => Kc,

        "BranchAndBoundClustering:max_abs_optimality_gap" => 0.,
        "BranchAndBoundClustering:max_rel_optimality_gap" => 0.,
        "BranchAndBoundClustering:E1_bound_in_rate_bound" => false,
        "BranchAndBoundClustering:store_fathomed_subtree_sizes" => true,
    ],
]

srand(725242)
network =
    setup_random_large_scale_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        num_streams=simulation_params["d"],
        geography_size=simulation_params["geography_size"],
        MS_serving_BS_distance=simulation_params["MS_serving_BS_distance"])

raw_results =
    simulate_assignment_convergence(network, simulation_params)

println("-- Saving results")
save("convergence.jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "raw_results", raw_results)
