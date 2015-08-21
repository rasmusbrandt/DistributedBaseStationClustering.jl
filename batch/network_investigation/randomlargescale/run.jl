#!/usr/bin/env julia

include("../../../src/IAClustering.jl")
using IAClustering, CoordinatedPrecoding, Compat
include("../plot_params.jl")

##########################################################################
# General settings
seed = 8367353
start_time = strftime("%Y%m%dT%H%M%S", time())

##########################################################################
# Initial simulation params
initial_simulation_params = [
    "simulation_name" => "initial",
    "I" => 8, "Kc" => 1, "N" => 2, "M" => 2, "d" => 1,
    "Ndrops" => 10, "Nsim" => 5,
    "geography_size" => (1500.,1500.),
    "MS_serving_BS_distance" => Nullable(150.),
    "assignment_methods" => [
        # ExhaustiveSearchClustering,
        # BranchAndBoundClustering,

        CoalitionFormationClustering_AttachOrSupplant,
        CoalitionFormationClustering_Attach,

        GreedyClustering_Single,
        GreedyClustering_Multiple,

        # Chen2014_LinearObj_ExhaustiveSearch,
        # Peters2012_Heuristic,

        GrandCoalitionClustering,
        RandomClustering,
        NoClustering,
    ],
    "precoding_methods" => [ RobustIntraclusterWMMSE, ],
    "aux_network_params" => [
        "num_coherence_symbols" => 2700,
        "beta_network_sdma" => 0.8,
    ],
    # "aux_assignment_params" => [ ],
    "aux_precoding_params" => [
        "initial_precoders" => "eigendirection",
        "stop_crit" => 1e-2,
        "max_iters" => 1000,
    ],
    "independent_variable" => (set_average_SNRs_dB!, -10:10:50),
]

##########################################################################
# Small scenario
simulation_params = deepcopy(initial_simulation_params)

network =
    setup_random_large_scale_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        num_streams=simulation_params["d"],
        geography_size=initial_simulation_params["geography_size"],
        MS_serving_BS_distance=initial_simulation_params["MS_serving_BS_distance"])

# Include high complexity methods for these simulations only
unshift!(simulation_params["assignment_methods"], BranchAndBoundClustering)

simulation_params["simulation_name"] = "small"

srand(seed)
raw_precoding_results, raw_assignment_results =
    simulate(network, simulation_params, loop_over=:assignment_methods)

processed_results = postprocess(raw_precoding_results, simulation_params, plot_params_instantaneous_sumrate)
plot(processed_results, simulation_params, plot_params_instantaneous_sumrate)
processed_results = postprocess(raw_assignment_results, simulation_params, plot_params_longterm_sumrate)
plot(processed_results, simulation_params, plot_params_longterm_sumrate)
processed_results = postprocess(raw_assignment_results, simulation_params, plot_params_longterm_avg_cluster_size)
plot(processed_results, simulation_params, plot_params_longterm_avg_cluster_size)
processed_results = postprocess(raw_assignment_results, simulation_params, plot_params_longterm_num_sum_throughput_calculations)
plot(processed_results, simulation_params, plot_params_longterm_num_sum_throughput_calculations)

##########################################################################
# Large scenario
simulation_params = deepcopy(initial_simulation_params)

simulation_params["I"] = 2*initial_simulation_params["I"]

network =
    setup_random_large_scale_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        num_streams=simulation_params["d"],
        geography_size=initial_simulation_params["geography_size"],
        MS_serving_BS_distance=initial_simulation_params["MS_serving_BS_distance"])

simulation_params["simulation_name"] = "large"

srand(seed)
raw_precoding_results, raw_assignment_results =
    simulate(network, simulation_params, loop_over=:assignment_methods)

processed_results = postprocess(raw_precoding_results, simulation_params, plot_params_instantaneous_sumrate)
plot(processed_results, simulation_params, plot_params_instantaneous_sumrate)
processed_results = postprocess(raw_assignment_results, simulation_params, plot_params_longterm_sumrate)
plot(processed_results, simulation_params, plot_params_longterm_sumrate)
processed_results = postprocess(raw_assignment_results, simulation_params, plot_params_longterm_avg_cluster_size)
plot(processed_results, simulation_params, plot_params_longterm_avg_cluster_size)
processed_results = postprocess(raw_assignment_results, simulation_params, plot_params_longterm_num_sum_throughput_calculations)
plot(processed_results, simulation_params, plot_params_longterm_num_sum_throughput_calculations)
