#!/usr/bin/env julia

require("../../../src/LongtermIAClustering.jl")
using LongtermIAClustering, CoordinatedPrecoding, Compat
require("../plot_params.jl")

##########################################################################
# General settings
seed = 7256141
start_time = strftime("%Y%m%dT%H%M%S", time())
SNRs_dB = -20:10:50

##########################################################################
# Initial simulation params
initial_simulation_params = [
    "simulation_name" => "initial",
    "Ic" => 2, "Kc" => 3, "N" => 2, "M" => 2, "d" => 1,
    "pico_centre_distance" => 100.,
    "Ndrops" => 10, "Nsim" => 5,
    "assignment_methods" => [
        # ExhaustiveSearchClustering,
        BranchAndBoundClustering,

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
    "aux_assignment_params" => [
        "max_num_MSs_per_BS" => 1,
    ],
    "aux_precoding_params" => [
        "initial_precoders" => "eigendirection",
        "stop_crit" => 1e-2,
        "max_iters" => 1000,
    ],
    "independent_variable" => (set_transmit_powers_dBm!, [ [p*ones(3), (p-20)*ones(6)] for p = SNRs_dB ]),
]

plot_params_instantaneous_sumrate["xvals"] = SNRs_dB
plot_params_longterm_sumrate["xvals"] = SNRs_dB
plot_params_longterm_avg_cluster_size["xvals"] = SNRs_dB

##########################################################################
# Small scenario
simulation_params = deepcopy(initial_simulation_params)

network =
    setup_triangularhetnet_network(simulation_params["Ic"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        num_streams=simulation_params["d"],
        pico_centre_distance=simulation_params["pico_centre_distance"])

simulation_params["simulation_name"] = ""

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
