#!/usr/bin/env julia

require("../../../../../src/IAClustering.jl")
using IAClustering, CoordinatedPrecoding
require("../../../plot_params.jl")

##########################################################################
# General settings
seed = 28373636
start_time = strftime("%Y%m%dT%H%M%S", time())

##########################################################################
# Initial simulation params
initial_simulation_params = [
    "simulation_name" => "initial",
    "I" => 8, "Kc" => 1, "N" => 2, "M" => 2, "d" => 1,
    "Ndrops" => 10, "Nsim" => 5,
    "geography_size" => (250.,250.),
    "MS_serving_BS_distance" => Nullable{Float64}(),
    "assignment_methods" => [
        # ExhaustiveSearchClustering,
        # BranchAndBoundClustering,

        CoalitionFormationClustering_Group,
        CoalitionFormationClustering_Individual,

        GreedyClustering_Single,
        GreedyClustering_Multiple,

        # Chen2014_LinearObj_ExhaustiveSearch,
        # Peters2012_Heuristic,

        GrandCoalitionClustering,
        RandomClustering,
        NoClustering,
    ],
    "precoding_methods" => [
        RobustIntraclusterWMMSE,
    ],
    "aux_network_params" => [
        "num_coherence_symbols" => 2160,
    ],
    "aux_assignment_params" => [
        "clustering_type" => :spectrum_sharing,

        "apply_overhead_prelog" => true,
        "IA_infeasible_negative_inf_utility" => true,
        "replace_E1_utility_with_lower_bound" => false,
    ],
    "aux_precoding_params" => [
        "initial_precoders" => "eigendirection",
        "stop_crit" => 1e-3,
        "max_iters" => 1000,
    ],
    "independent_variable" => (set_transmit_powers_dBm!, -80:10:0),
]

##########################################################################
# Small scenario with/without overhead
simulation_params = deepcopy(initial_simulation_params)

network =
    setup_random_large_scale_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        num_streams=simulation_params["d"],
        geography_size=simulation_params["geography_size"],
        MS_serving_BS_distance=simulation_params["MS_serving_BS_distance"])

# Include high complexity methods for these simulations only
unshift!(simulation_params["assignment_methods"], Chen2014_LinearObj_ExhaustiveSearch)
unshift!(simulation_params["assignment_methods"], BranchAndBoundClustering)
unshift!(simulation_params["assignment_methods"], ExhaustiveSearchClustering)

simulation_params["simulation_name"] = "small-with_overhead"
simulation_params["aux_assignment_params"]["apply_overhead_prelog"] = true

srand(seed)
raw_precoding_results, raw_assignment_results =
    simulate(network, simulation_params, loop_over=:assignment_methods)

processed_results = postprocess(raw_precoding_results, simulation_params, plot_params_instantaneous_sumrate)
plot(processed_results, simulation_params, plot_params_instantaneous_sumrate)
processed_results = postprocess(raw_assignment_results, simulation_params, plot_params_longterm_sumrate)
plot(processed_results, simulation_params, plot_params_longterm_sumrate)
processed_results = postprocess(raw_assignment_results, simulation_params, plot_params_longterm_num_utility_calculations)
plot(processed_results, simulation_params, plot_params_longterm_num_utility_calculations)
processed_results = postprocess(raw_assignment_results, simulation_params, plot_params_longterm_num_longterm_rate_calculations)
plot(processed_results, simulation_params, plot_params_longterm_num_longterm_rate_calculations)
processed_results = postprocess(raw_assignment_results, simulation_params, plot_params_longterm_clusters)
plot(processed_results, simulation_params, plot_params_longterm_clusters)

simulation_params["simulation_name"] = "small-without_overhead"
simulation_params["aux_assignment_params"]["apply_overhead_prelog"] = false

srand(seed)
raw_precoding_results, raw_assignment_results =
    simulate(network, simulation_params, loop_over=:assignment_methods)

processed_results = postprocess(raw_precoding_results, simulation_params, plot_params_instantaneous_sumrate)
plot(processed_results, simulation_params, plot_params_instantaneous_sumrate)
processed_results = postprocess(raw_assignment_results, simulation_params, plot_params_longterm_sumrate)
plot(processed_results, simulation_params, plot_params_longterm_sumrate)
processed_results = postprocess(raw_assignment_results, simulation_params, plot_params_longterm_num_utility_calculations)
plot(processed_results, simulation_params, plot_params_longterm_num_utility_calculations)
processed_results = postprocess(raw_assignment_results, simulation_params, plot_params_longterm_num_longterm_rate_calculations)
plot(processed_results, simulation_params, plot_params_longterm_num_longterm_rate_calculations)
processed_results = postprocess(raw_assignment_results, simulation_params, plot_params_longterm_clusters)
plot(processed_results, simulation_params, plot_params_longterm_clusters)

##########################################################################
# Large scenario with/without overhead
simulation_params = deepcopy(initial_simulation_params)

simulation_params["I"] = 2*initial_simulation_params["I"]

network =
    setup_random_large_scale_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        num_streams=simulation_params["d"],
        geography_size=simulation_params["geography_size"])

simulation_params["simulation_name"] = "large-with_overhead"
simulation_params["aux_assignment_params"]["apply_overhead_prelog"] = true

srand(seed)
raw_precoding_results, raw_assignment_results =
    simulate(network, simulation_params, loop_over=:assignment_methods)

processed_results = postprocess(raw_precoding_results, simulation_params, plot_params_instantaneous_sumrate)
plot(processed_results, simulation_params, plot_params_instantaneous_sumrate)
processed_results = postprocess(raw_assignment_results, simulation_params, plot_params_longterm_sumrate)
plot(processed_results, simulation_params, plot_params_longterm_sumrate)
processed_results = postprocess(raw_assignment_results, simulation_params, plot_params_longterm_num_utility_calculations)
plot(processed_results, simulation_params, plot_params_longterm_num_utility_calculations)
processed_results = postprocess(raw_assignment_results, simulation_params, plot_params_longterm_num_longterm_rate_calculations)
plot(processed_results, simulation_params, plot_params_longterm_num_longterm_rate_calculations)
processed_results = postprocess(raw_assignment_results, simulation_params, plot_params_longterm_clusters)
plot(processed_results, simulation_params, plot_params_longterm_clusters)

simulation_params["simulation_name"] = "large-without_overhead"
simulation_params["aux_assignment_params"]["apply_overhead_prelog"] = false

srand(seed)
raw_precoding_results, raw_assignment_results =
    simulate(network, simulation_params, loop_over=:assignment_methods)

processed_results = postprocess(raw_precoding_results, simulation_params, plot_params_instantaneous_sumrate)
plot(processed_results, simulation_params, plot_params_instantaneous_sumrate)
processed_results = postprocess(raw_assignment_results, simulation_params, plot_params_longterm_sumrate)
plot(processed_results, simulation_params, plot_params_longterm_sumrate)
processed_results = postprocess(raw_assignment_results, simulation_params, plot_params_longterm_num_utility_calculations)
plot(processed_results, simulation_params, plot_params_longterm_num_utility_calculations)
processed_results = postprocess(raw_assignment_results, simulation_params, plot_params_longterm_num_longterm_rate_calculations)
plot(processed_results, simulation_params, plot_params_longterm_num_longterm_rate_calculations)
processed_results = postprocess(raw_assignment_results, simulation_params, plot_params_longterm_clusters)
plot(processed_results, simulation_params, plot_params_longterm_clusters)
