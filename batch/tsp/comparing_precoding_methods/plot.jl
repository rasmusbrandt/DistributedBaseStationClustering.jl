#!/usr/bin/env julia

include(joinpath(dirname(@__FILE__), "../../../src/DistributedBSClustering.jl"))
using CoordinatedPrecoding, DistributedBSClustering
using Compat, JLD, LaTeXStrings

include(joinpath(dirname(@__FILE__), "../simulation_params.jl"))
include(joinpath(dirname(@__FILE__), "../simulation_params-precoding_methods.jl"))
include(joinpath(dirname(@__FILE__), "../simulation_params-SNR.jl"))
include(joinpath(dirname(@__FILE__), "../plot_params-final.jl"))

##########################################################################
# Setup
sim_name = "precoding_methods"
data = load("$(sim_name).jld")

# Defaults
PyPlot.rc("lines", linewidth=1., markersize=3, markeredgewidth=0.5)
PyPlot.rc("font", size=8, family="serif", serif="Computer Modern Sans Serif")
PyPlot.rc("text", usetex=true)
PyPlot.rc("text.latex", preamble="\\usepackage{amsmath}")
PyPlot.rc("axes", linewidth=0.5, labelsize=8)
PyPlot.rc("xtick", labelsize=8)
PyPlot.rc("ytick", labelsize=8)
PyPlot.rc("legend", fancybox=true, fontsize=6)
PyPlot.rc("figure", figsize=(3.50,2.16), dpi=125)

# Legend creation helper
show_legend!(ax, loc="lower center") = begin
    legend = ax[:legend](loc=loc)
    legend_frame = legend[:get_frame]()
    PyPlot.setp(legend_frame, linewidth=0.5)
end

# Extract the values easily
xvals = data["simulation_params"]["independent_variable"][2]
yvals_precoding(method, plot_val) = data["processed_precoding_results"][mean_idx][method][plot_val]

##########################################################################
# instantaneous-sumrate
plot_name = "instantaneous-sumrate"
fig = PyPlot.figure()
ax = fig[:add_axes]((0.12,0.15,0.98-0.12,0.96-0.15))
ax[:plot](xvals, yvals_precoding("RobustIntraclusterWMMSE", "weighted_logdet_rates_full"),
    color=colours_precoding[:RobustIntraclusterWMMSE],
    linestyle="-",
    marker=markers_precoding[:RobustIntraclusterWMMSE], markeredgecolor=colours_precoding[:RobustIntraclusterWMMSE], markevery=5,
    label="Robust WMMSE (ICI aware)")
ax[:plot](xvals, yvals_precoding("RobustIntraclusterWMMSE", "weighted_logdet_rates_partial"),
    color=colours_precoding[:RobustIntraclusterWMMSE],
    linestyle="--",
    marker=markers_precoding[:RobustIntraclusterWMMSE], markeredgecolor=colours_precoding[:RobustIntraclusterWMMSE], markevery=5,
    label="Robust WMMSE (ICI oblivious)")
ax[:plot](xvals, yvals_precoding("RobustChen2014_MaxSINR", "weighted_logdet_rates_full"),
    color=colours_precoding[:RobustChen2014_MaxSINR],
    linestyle="-",
    marker=markers_precoding[:RobustChen2014_MaxSINR], markeredgecolor=colours_precoding[:RobustChen2014_MaxSINR], markevery=5,
    label="Robust MaxSINR [x] (ICI aware)")
ax[:plot](xvals, yvals_precoding("RobustChen2014_MaxSINR", "weighted_logdet_rates_partial"),
    color=colours_precoding[:RobustChen2014_MaxSINR],
    linestyle="--",
    marker=markers_precoding[:RobustChen2014_MaxSINR], markeredgecolor=colours_precoding[:RobustChen2014_MaxSINR], markevery=5,
    label="Robust MaxSINR [x] (ICI oblivious)")
ax[:plot](xvals, yvals_precoding("NaiveIntraclusterWMMSE", "weighted_logdet_rates_full"),
    color=colours_precoding[:NaiveIntraclusterWMMSE],
    linestyle="-",
    marker=markers_precoding[:NaiveIntraclusterWMMSE], markeredgecolor=colours_precoding[:NaiveIntraclusterWMMSE], markevery=5,
    label="Naive WMMSE (ICI aware)")
ax[:plot](xvals, yvals_precoding("NaiveIntraclusterWMMSE", "weighted_logdet_rates_partial"),
    color=colours_precoding[:NaiveIntraclusterWMMSE],
    linestyle="--",
    marker=markers_precoding[:NaiveIntraclusterWMMSE], markeredgecolor=colours_precoding[:NaiveIntraclusterWMMSE], markevery=5,
    label="Naive WMMSE (ICI oblivious)")

ax[:set_ylim]([-2, 100])
ax[:set_xlabel]("Signal-to-noise ratio [dB]")
ax[:set_ylabel]("Sum throughput [bits/s/Hz]")
show_legend!(ax, "lower right")
fig[:savefig]("$(sim_name)_$(plot_name).eps")
