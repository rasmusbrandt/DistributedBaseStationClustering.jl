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
data = load("raw-large_network.jld")

##########################################################################
# Perform post processing
postprocess_params_assignment = [
    "objective" => :sumrate,
    "methods" => [
        "CoalitionFormationClustering_Individual" => [
            ("utilities",),
        ],

        "GrandCoalitionClustering" => [
            ("utilities",),
        ],

        "RandomClustering" => [
            ("utilities",),
        ],

        "NoClustering" => [
            ("utilities",),
        ],
    ]
]
results_assignment, results_assignment_mean, results_assignment_var = postprocess(data["raw_assignment_results"], data["simulation_params"], postprocess_params_assignment)

postprocess_params_precoding = [
    "objective" => :sumrate,
    "methods" => [
        "CoalitionFormationClustering_Individual" => [
            ("weighted_logdet_rates_LB",),
        ],

        "GrandCoalitionClustering" => [
            ("weighted_logdet_rates_LB",),
        ],

        "RandomClustering" => [
            ("weighted_logdet_rates_LB",),
        ],

        "NoClustering" => [
            ("weighted_logdet_rates_LB",),
        ],
    ]
]
results_precoding, results_precoding_mean, results_precoding_var = postprocess(data["raw_precoding_results"], data["simulation_params"], postprocess_params_precoding)

##########################################################################
# Build figures
PyPlot.rc("lines", linewidth=1, markersize=6)
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

ax1[:plot](transmit_powers_dBm, results_assignment_mean["CoalitionFormationClustering_Individual"]["utilities"][:,2], color=cf_col, linestyle="-", marker=".", label=L"Coalition formation ($b_k = 100$)")
ax1[:plot](transmit_powers_dBm, results_assignment_mean["CoalitionFormationClustering_Individual"]["utilities"][:,1], color=cf_col, linestyle="--", marker=".", label=L"Coalition formation ($b_k = 10$)")
ax1[:plot](transmit_powers_dBm, results_assignment_mean["RandomClustering"]["utilities"][:,1], color=random_col, linestyle="-", marker=".", label="Random coalitions")
ax1[:plot](transmit_powers_dBm, results_assignment_mean["NoClustering"]["utilities"][:,1], color=singleton_col, linestyle="-", marker=".", label="Singleton coalitions")
ax1[:plot](transmit_powers_dBm, zeros(results_assignment_mean["GrandCoalitionClustering"]["utilities"][:,1]), color=grand_col, linestyle="-", marker=".", label="Grand coalition")

ax1[:set_ylim](-5, 65)

ax1[:set_xlabel]("Transmit power [dBm]")
ax1[:set_ylabel]("Longterm sum throughput [bits/s/Hz]")

legend = ax1[:legend](loc="lower right")
# legend_lines = legend[:get_lines]()
legend_frame = legend[:get_frame]()
# PyPlot.setp(legend_lines, linewidth=0.5)
PyPlot.setp(legend_frame, linewidth=0.5)


fig3 = PyPlot.figure()
ax3 = fig3[:add_axes]((0.11,0.16,0.95-0.11,0.95-0.16))

ax3[:plot](transmit_powers_dBm, results_precoding_mean["CoalitionFormationClustering_Individual"]["weighted_logdet_rates_LB"][:,2], color=cf_col, linestyle="-", marker=".", label=L"Coalition formation ($b_k = 100$)")
ax3[:plot](transmit_powers_dBm, results_precoding_mean["CoalitionFormationClustering_Individual"]["weighted_logdet_rates_LB"][:,1], color=cf_col, linestyle="--", marker=".", label=L"Coalition formation ($b_k = 10$)")
ax3[:plot](transmit_powers_dBm, results_precoding_mean["RandomClustering"]["weighted_logdet_rates_LB"][:,1], color=random_col, linestyle="-", marker=".", label="Random coalitions")
ax3[:plot](transmit_powers_dBm, results_precoding_mean["NoClustering"]["weighted_logdet_rates_LB"][:,1], color=singleton_col, linestyle="-", marker=".", label="Singleton coalitions")
ax3[:plot](transmit_powers_dBm, results_precoding_mean["GrandCoalitionClustering"]["weighted_logdet_rates_LB"][:,1], color=grand_col, linestyle="-", marker=".", label="Grand coalition")

ax3[:set_ylim](-5, 65)

ax3[:set_xlabel]("Transmit power [dBm]")
ax3[:set_ylabel]("Avg. instant. sum rate [bits/s/Hz]")

legend = ax3[:legend](loc="best")
# legend_lines = legend[:get_lines]()
legend_frame = legend[:get_frame]()
# PyPlot.setp(legend_lines, linewidth=0.5)
PyPlot.setp(legend_frame, linewidth=0.5)


##########################################################################
# Write files
fig1[:savefig]("large_network-longterm-sumrate.pdf")
fig3[:savefig]("large_network-instantaneous-sumrate.pdf")
