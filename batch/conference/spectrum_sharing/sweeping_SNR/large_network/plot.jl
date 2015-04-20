#!/usr/bin/env julia

include(joinpath(dirname(@__FILE__), "../../../../../src/IAClustering.jl"))
using IAClustering, CoordinatedPrecoding
using HDF5, JLD
using LaTeXStrings

include(joinpath(dirname(@__FILE__), "../../../simulation_params.jl"))
include(joinpath(dirname(@__FILE__), "../../../simulation_params-assignment_methods.jl"))
include(joinpath(dirname(@__FILE__), "../../../simulation_params-large_network.jl"))
include(joinpath(dirname(@__FILE__), "../../../simulation_params-SNR.jl"))
include(joinpath(dirname(@__FILE__), "../../../plot_params-assignment_methods.jl"))

##########################################################################
# Load data
using HDF5, JLD
data = load("raw-large_network-WMMSE.jld")

##########################################################################
# Perform post processing
postprocess_params_assignment = [
    "objective" => :sumrate,
    "methods" => [
        "CoalitionFormationClustering_Individual" => [
            ("utilities",),
            ("no_clusters",),
        ],

        "GreedyClustering_Multiple" => [
            ("utilities",),
            ("no_clusters",),
        ],

        "RandomClustering" => [
            ("utilities",),
            ("no_clusters",),
        ],

        "NoClustering" => [
            ("utilities",),
            ("no_clusters",),
        ],
    ]
]
results_assignment, results_assignment_mean, results_assignment_var = postprocess(data["raw_assignment_results"], data["simulation_params"], postprocess_params_assignment)

postprocess_params_precoding = [
    "objective" => :sumrate,
    "methods" => [
        "CoalitionFormationClustering_Individual" => [
            ("weighted_logdet_rates_partial",),
        ],

        "GreedyClustering_Multiple" => [
            ("weighted_logdet_rates_partial",),
        ],

        "GrandCoalitionClustering" => [
            ("weighted_logdet_rates_partial",),
        ],

        "RandomClustering" => [
            ("weighted_logdet_rates_partial",),
        ],

        "NoClustering" => [
            ("weighted_logdet_rates_partial",),
        ],
    ]
]
results_precoding, results_precoding_mean, results_precoding_var = postprocess(data["raw_precoding_results"], data["simulation_params"], postprocess_params_precoding)

##########################################################################
# Build figures
PyPlot.rc("lines", linewidth=1., markersize=6)
PyPlot.rc("font", size=8, family="serif", serif="Computer Modern Sans Serif")
PyPlot.rc("text", usetex=true)
PyPlot.rc("text.latex", preamble="\\usepackage{amsmath}")
PyPlot.rc("pgf", preamble="\\usepackage{amsmath}")
PyPlot.rc("axes", linewidth=0.5, labelsize=8)
PyPlot.rc("xtick", labelsize=8)
PyPlot.rc("ytick", labelsize=8)
PyPlot.rc("legend", fancybox=true, fontsize=6)
PyPlot.rc("figure", figsize=(3.40,2.16), dpi=125)

fig1 = PyPlot.figure()
ax1 = fig1[:add_axes]((0.11,0.15,0.95-0.11,0.95-0.15))

ax1[:plot](transmit_powers_dBm, results_assignment_mean["CoalitionFormationClustering_Individual"]["utilities"], color="LimeGreen", linestyle="-", label="Distributed coalition formation")
ax1[:plot](transmit_powers_dBm, results_assignment_mean["GreedyClustering_Multiple"]["utilities"], color="DarkOrchid", linestyle="-", label="Centralized greedy algorithm")
ax1[:plot](transmit_powers_dBm, results_assignment_mean["RandomClustering"]["utilities"], color="Khaki", linestyle="-", label="Random IA feasible coalitions")
ax1[:plot](transmit_powers_dBm, results_assignment_mean["NoClustering"]["utilities"], color="Pink", linestyle="-", label="Singleton coalitions")

ax1[:set_ylim](0, 60)

ax1[:set_xlabel]("Transmit power [dBm]")
ax1[:set_ylabel]("Longterm sum rate [bits/s/Hz]")

legend = ax1[:legend](loc="lower right")
# legend_lines = legend[:get_lines]()
legend_frame = legend[:get_frame]()
# PyPlot.setp(legend_lines, linewidth=0.5)
PyPlot.setp(legend_frame, linewidth=0.5)


fig2 = PyPlot.figure()
ax2 = fig2[:add_axes]((0.11,0.15,0.95-0.11,0.95-0.15))

ax2[:plot](transmit_powers_dBm, results_assignment_mean["CoalitionFormationClustering_Individual"]["no_clusters"], color="LimeGreen", linestyle="-", label="Distributed coalition formation")
ax2[:plot](transmit_powers_dBm, results_assignment_mean["GreedyClustering_Multiple"]["no_clusters"], color="DarkOrchid", linestyle="-", label="Centralized greedy algorithm")
ax2[:plot](transmit_powers_dBm, results_assignment_mean["RandomClustering"]["no_clusters"], color="Khaki", linestyle="-", label="Random IA feasible coalitions")
ax2[:plot](transmit_powers_dBm, results_assignment_mean["NoClustering"]["no_clusters"], color="Pink", linestyle="-", label="Singleton coalitions")

ax2[:set_ylim](0, 30)

ax2[:set_xlabel]("Transmit power [dBm]")
ax2[:set_ylabel]("Average number of coalitions formed")

legend = ax2[:legend](loc="best")
# legend_lines = legend[:get_lines]()
legend_frame = legend[:get_frame]()
# PyPlot.setp(legend_lines, linewidth=0.5)
PyPlot.setp(legend_frame, linewidth=0.5)


fig3 = PyPlot.figure()
ax3 = fig3[:add_axes]((0.11,0.15,0.95-0.11,0.95-0.15))

ax3[:plot](transmit_powers_dBm, results_precoding_mean["CoalitionFormationClustering_Individual"]["weighted_logdet_rates_partial"], color="LimeGreen", linestyle="-", label="Distributed coalition formation")
ax3[:plot](transmit_powers_dBm, results_precoding_mean["GreedyClustering_Multiple"]["weighted_logdet_rates_partial"], color="DarkOrchid", linestyle="-", label="Centralized greedy algorithm")
ax3[:plot](transmit_powers_dBm, results_precoding_mean["RandomClustering"]["weighted_logdet_rates_partial"], color="Khaki", linestyle="-", label="Random IA feasible coalitions")
ax3[:plot](transmit_powers_dBm, results_precoding_mean["NoClustering"]["weighted_logdet_rates_partial"], color="Pink", linestyle="-", label="Singleton coalitions")
ax3[:plot](transmit_powers_dBm, results_precoding_mean["GrandCoalitionClustering"]["weighted_logdet_rates_partial"], color="Maroon", linestyle="-", label="Grand coalition")

ax3[:set_ylim](-10, 100)

ax3[:set_xlabel]("Transmit power [dBm]")
ax3[:set_ylabel]("Average sum rate [bits/s/Hz]")

legend = ax3[:legend](loc="best")
# legend_lines = legend[:get_lines]()
legend_frame = legend[:get_frame]()
# PyPlot.setp(legend_lines, linewidth=0.5)
PyPlot.setp(legend_frame, linewidth=0.5)


##########################################################################
# Write files
fig1[:savefig]("large_network-longterm-sumrate.pgf")
fig1[:savefig]("large_network-longterm-sumrate.pdf")
fig2[:savefig]("large_network-longterm-no_clusters.pgf")
fig2[:savefig]("large_network-longterm-no_clusters.pdf")
fig3[:savefig]("large_network-instantaneous-sumrate.pgf")
fig3[:savefig]("large_network-instantaneous-sumrate.pdf")
