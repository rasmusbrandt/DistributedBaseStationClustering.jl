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
data = load("raw-small_network.jld")

##########################################################################
# Perform post processing
postprocess_params = [
    "objective" => :sumrate,
    "methods" => [
        "ExhaustiveSearchClustering" => [
            ("utilities",),
            ("no_clusters",),
        ],

        "Chen2014_ExhaustiveSearch" => [
            ("utilities",),
            ("no_clusters",),
        ],

        "CoalitionFormationClustering_Individual" => [
            ("utilities",),
            ("no_clusters",),
        ],

        "GrandCoalitionClustering" => [
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
results, results_mean, results_var = postprocess(data["raw_assignment_results"], data["simulation_params"], postprocess_params)

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

ax1[:plot](transmit_powers_dBm, results_mean["ExhaustiveSearchClustering"]["utilities"], color="Coral", linestyle="-", label="Exhaustive search (proposed utility model)")
ax1[:plot](transmit_powers_dBm, results_mean["Chen2014_ExhaustiveSearch"]["utilities"], color="DodgerBlue", linestyle="-", label="Exhaustive search (utility model in [Chen2014])")
ax1[:plot](transmit_powers_dBm, results_mean["CoalitionFormationClustering_Individual"]["utilities"], color="LimeGreen", linestyle="-", label="Distributed coalition formation")
ax1[:plot](transmit_powers_dBm, results_mean["RandomClustering"]["utilities"], color="Khaki", linestyle="-", label="Random IA feasible coalitions")
ax1[:plot](transmit_powers_dBm, results_mean["NoClustering"]["utilities"], color="Pink", linestyle="-", label="Singleton coalitions")

ax1[:set_ylim](0, 25)

ax1[:set_xlabel]("Transmit power [dBm]")
ax1[:set_ylabel]("Longterm sum rate [bits/s/Hz]")

legend = ax1[:legend](loc="upper left")
# legend_lines = legend[:get_lines]()
legend_frame = legend[:get_frame]()
# PyPlot.setp(legend_lines, linewidth=0.5)
PyPlot.setp(legend_frame, linewidth=0.5)


fig2 = PyPlot.figure()
ax2 = fig2[:add_axes]((0.11,0.15,0.95-0.11,0.95-0.15))

ax2[:plot](transmit_powers_dBm, results_mean["ExhaustiveSearchClustering"]["no_clusters"], color="Coral", linestyle="-", label="Exhaustive search (proposed utility model)")
ax2[:plot](transmit_powers_dBm, results_mean["Chen2014_ExhaustiveSearch"]["no_clusters"], color="DodgerBlue", linestyle="-", label="Exhaustive search (utility model in [Chen2014])")
ax2[:plot](transmit_powers_dBm, results_mean["CoalitionFormationClustering_Individual"]["no_clusters"], color="LimeGreen", linestyle="-", label="Distributed coalition formation")
ax2[:plot](transmit_powers_dBm, results_mean["RandomClustering"]["no_clusters"], color="Khaki", linestyle="-", label="Random IA feasible coalitions")
ax2[:plot](transmit_powers_dBm, results_mean["NoClustering"]["no_clusters"], color="Pink", linestyle="-", label="Singleton coalitions")

ax2[:set_ylim](0, 9)

ax2[:set_xlabel]("Transmit power [dBm]")
ax2[:set_ylabel]("Average number of coalitions formed")

legend = ax2[:legend](loc="upper right")
# legend_lines = legend[:get_lines]()
legend_frame = legend[:get_frame]()
# PyPlot.setp(legend_lines, linewidth=0.5)
PyPlot.setp(legend_frame, linewidth=0.5)

##########################################################################
# Write files
fig1[:savefig]("small_network-longterm-sumrate.pgf")
fig1[:savefig]("small_network-longterm-sumrate.pdf")
fig2[:savefig]("small_network-longterm-no_clusters.pgf")
fig2[:savefig]("small_network-longterm-no_clusters.pdf")
