#!/usr/bin/env julia

include(joinpath(dirname(@__FILE__), "../../../src/IAClustering.jl"))
using CoordinatedPrecoding, IAClustering
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
plot_name = "instantaneous-sumrate"; plot_val = "weighted_logdet_rates_full"
fig = PyPlot.figure()
ax = fig[:add_axes]((0.11,0.15,0.95-0.11,0.95-0.15))

for method in [:RobustIntraclusterWMMSE, :NaiveIntraclusterWMMSE, :RobustChen2014_MaxSINR, :NaiveChen2014_MaxSINR]
    ax[:plot](xvals, yvals_precoding(string(method), "weighted_logdet_rates_full")[:,2],
        color=colours_precoding[method],
        linestyle="-",
        marker=markers_precoding[method], markeredgecolor=colours_precoding[method], markevery=5,
        label=labels_precoding[method])
    ax[:plot](xvals, yvals_precoding(string(method), "weighted_logdet_rates_partial")[:,2],
        color=colours_precoding[method],
        linestyle="--",
        marker=markers_precoding[method], markeredgecolor=colours_precoding[method], markevery=5,
        label=labels_precoding[method])
end
ax[:set_ylim]([-2, 100])
ax[:set_xlabel]("Signal-to-noise ratio [dB]")
ax[:set_ylabel]("Average sum throughput [bits/s/Hz]")
show_legend!(ax, "upper left")
fig[:savefig]("$(sim_name)_$(plot_name).eps")
