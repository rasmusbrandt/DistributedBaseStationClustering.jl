#!/usr/bin/env julia

include("../../../../src/IAClustering.jl")
using IAClustering, CoordinatedPrecoding
include("../plot_params.jl")

##########################################################################
# General settings
srand(725242)
start_time = strftime("%Y%m%dT%H%M%S", time())
SNRs_dB = -20:10:50

##########################################################################
# Initial simulation params
initial_simulation_params = [
    "simulation_name" => "initial",
    "Ic" => 2, "Kc" => 3, "N" => 2, "M" => 2, "d" => 1,
    "pico_centre_distance" => 150.,
    "Ndrops" => 10, "Nsim" => 10,
    "assignment_methods" => [
        # ExhaustiveSearchClustering,
        # BranchAndBoundClustering,

        CoalitionFormationClustering_Group,
        CoalitionFormationClustering_Individual,

        # Chen2014_ExhaustiveSearch,
        # Peters2012_Heuristic,

        GrandCoalitionClustering,
        GreedyClustering,
        RandomClustering,
        NoClustering,
    ],
    "precoding_methods" => [
        RobustIntraclusterWMMSE,
    ],
    "aux_network_params" => [
        "no_coherence_symbols" => 1000,
    ],
    "aux_assignment_params" => [
        "max_MSs_per_BS" => 1, 
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
    "independent_variable" => (set_transmit_powers_dBm!, [ [p*ones(3), (p-20)*ones(6)] for p = SNRs_dB ]),
]

plot_params_instantaneous_sumrate["xvals"] = SNRs_dB
plot_params_longterm_sumrate["xvals"] = SNRs_dB
plot_params_longterm_iters["xvals"] = SNRs_dB

##########################################################################
# Small scenario with/without overhead
simulation_params = deepcopy(initial_simulation_params)

network =
    setup_triangularhetnet_network(simulation_params["Ic"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        no_streams=simulation_params["d"],
        pico_centre_distance=simulation_params["pico_centre_distance"])

# Include high complexity methods for these simulations only
unshift!(simulation_params["assignment_methods"], BranchAndBoundClustering)
unshift!(simulation_params["assignment_methods"], ExhaustiveSearchClustering)

simulation_params["simulation_name"] = "small-with_overhead"
simulation_params["aux_assignment_params"]["apply_overhead_prelog"] = true

raw_precoding_results, raw_assignment_results =
    simulate(network, simulation_params, loop_over=:assignment_methods)

processed_results = postprocess(raw_precoding_results, simulation_params, plot_params_instantaneous_sumrate)
plot(processed_results, simulation_params, plot_params_instantaneous_sumrate)
processed_results = postprocess(raw_assignment_results, simulation_params, plot_params_longterm_sumrate)
plot(processed_results, simulation_params, plot_params_longterm_sumrate)
processed_results = postprocess(raw_assignment_results, simulation_params, plot_params_longterm_iters)
plot(processed_results, simulation_params, plot_params_longterm_iters)

simulation_params["simulation_name"] = "small-without_overhead"
simulation_params["aux_assignment_params"]["apply_overhead_prelog"] = false

raw_precoding_results, raw_assignment_results =
    simulate(network, simulation_params, loop_over=:assignment_methods)

processed_results = postprocess(raw_precoding_results, simulation_params, plot_params_instantaneous_sumrate)
plot(processed_results, simulation_params, plot_params_instantaneous_sumrate)
processed_results = postprocess(raw_assignment_results, simulation_params, plot_params_longterm_sumrate)
plot(processed_results, simulation_params, plot_params_longterm_sumrate)
processed_results = postprocess(raw_assignment_results, simulation_params, plot_params_longterm_iters)
plot(processed_results, simulation_params, plot_params_longterm_iters)

