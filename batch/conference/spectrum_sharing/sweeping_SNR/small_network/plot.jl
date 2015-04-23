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
PyPlot.rc("lines", linewidth=1, markersize=5)
PyPlot.rc("font", size=8, family="serif", serif="Computer Modern Sans Serif")
PyPlot.rc("text", usetex=true)
PyPlot.rc("text.latex", preamble="\\usepackage{amsmath}")
PyPlot.rc("axes", linewidth=0.5, labelsize=8)
PyPlot.rc("xtick", labelsize=8)
PyPlot.rc("ytick", labelsize=8)
PyPlot.rc("legend", fancybox=true, fontsize=6)
PyPlot.rc("figure", figsize=(3.5,2.0))

optimal_col = "#e7298a"
cf_col = "#1b9e77"
chen2014_col = "#66a61e"
random_col = "#d95f02"
singleton_col = "#7570b3"
grand_col = "#e6ab02"

fig1 = PyPlot.figure()
ax1 = fig1[:add_axes]((0.11,0.16,0.95-0.11,0.95-0.16))

ax1[:plot](transmit_powers_dBm, results_mean["ExhaustiveSearchClustering"]["utilities"][:,1], color=optimal_col, linestyle="-", marker=".", label="Exhaustive search")
ax1[:plot](transmit_powers_dBm, results_mean["CoalitionFormationClustering_Individual"]["utilities"][:,1], color=cf_col, linestyle="-", marker=".", label=L"Coalition formation ($b_k = 10$)")
ax1[:plot](transmit_powers_dBm, results_mean["RandomClustering"]["utilities"][:,1], color=random_col, linestyle="-", marker=".", label="Random coalitions")
ax1[:plot](transmit_powers_dBm, results_mean["NoClustering"]["utilities"][:,1], color=singleton_col, linestyle="-", marker=".", label="Singleton coalitions")
ax1[:plot](transmit_powers_dBm, zeros(results_mean["GrandCoalitionClustering"]["utilities"][:,1]), color=grand_col, linestyle="-", marker=".", label="Grand coalition")

ax1[:set_ylim](-2, 25)

ax1[:set_xlabel]("Transmit power [dBm]")
ax1[:set_ylabel]("Longterm throughput [bits/s/Hz]")

legend = ax1[:legend](loc="lower right")
# legend_lines = legend[:get_lines]()
legend_frame = legend[:get_frame]()
# PyPlot.setp(legend_lines, linewidth=0.5)
PyPlot.setp(legend_frame, linewidth=0.5)


fig2 = PyPlot.figure()
ax2 = fig2[:add_axes]((0.11,0.16,0.95-0.11,0.95-0.16))

ax2[:plot](transmit_powers_dBm, data["simulation_params"]["I"]./results_mean["CoalitionFormationClustering_Individual"]["no_clusters"][:,1], color=cf_col, linestyle="-", marker=".", label=L"Coalition formation ($b_k = 10$)")
ax2[:plot](transmit_powers_dBm, data["simulation_params"]["I"]./results_mean["ExhaustiveSearchClustering"]["no_clusters"][:,1], color=optimal_col, linestyle="-", marker=".", label="Exhaustive search")
ax2[:plot](transmit_powers_dBm, data["simulation_params"]["I"]./results_mean["RandomClustering"]["no_clusters"][:,1], color=random_col, linestyle="-", marker=".", label="Random coalitions")
ax2[:plot](transmit_powers_dBm, data["simulation_params"]["I"]./results_mean["NoClustering"]["no_clusters"][:,1], color=singleton_col, linestyle="-", marker=".", label="Singleton coalitions")

ax2[:set_ylim](0.5, 3.2)
ax2[:set_yticks]([1, 2, 3])

ax2[:set_xlabel]("Transmit power [dBm]")
ax2[:set_ylabel]("Average size of coalitions")

legend = ax2[:legend](loc="lower right")
# legend_lines = legend[:get_lines]()
legend_frame = legend[:get_frame]()
# PyPlot.setp(legend_lines, linewidth=0.5)
PyPlot.setp(legend_frame, linewidth=0.5)

##########################################################################
# Write files
fig1[:savefig]("small_network-longterm-sumrate.pdf")
fig2[:savefig]("small_network-longterm-no_clusters.pdf")
