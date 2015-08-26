#!/usr/bin/env julia

include(joinpath(dirname(@__FILE__), "../../../src/IAClustering.jl"))
using CoordinatedPrecoding, IAClustering
using Compat, JLD, LaTeXStrings

include(joinpath(dirname(@__FILE__), "../simulation_params.jl"))
include(joinpath(dirname(@__FILE__), "../simulation_params-assignment_methods.jl"))
include(joinpath(dirname(@__FILE__), "../simulation_params-beta.jl"))
include(joinpath(dirname(@__FILE__), "../plot_params-final.jl"))

##########################################################################
# Setup
sim_name = "beta"
data = load("$(sim_name).jld")

# Defaults
PyPlot.rc("lines", linewidth=1., markersize=3, markeredgewidth=0.5)
PyPlot.rc("font", size=8, family="serif", serif="Computer Modern Sans Serif")
PyPlot.rc("text", usetex=true)
PyPlot.rc("text.latex", preamble="\\usepackage{amsmath}")
PyPlot.rc("axes", linewidth=0.5, labelsize=8)
PyPlot.rc("xtick", labelsize=8)
PyPlot.rc("ytick", labelsize=8)
PyPlot.rc("legend", fancybox=true, fontsize=7)
PyPlot.rc("figure", figsize=(3.50,2.16), dpi=125)

# Legend creation helper
show_legend!(ax, loc="lower center") = begin
    legend = ax[:legend](loc=loc)
    legend_frame = legend[:get_frame]()
    PyPlot.setp(legend_frame, linewidth=0.5)
end

# Extract the values easily
xvals = data["simulation_params"]["independent_variable"][2]
yvals_assignment(method, plot_val) = data["processed_assignment_results"][mean_idx][method][plot_val]
yvals_precoding(method, plot_val) = data["processed_precoding_results"][mean_idx][method][plot_val]

##########################################################################
# longterm-sumrate
plot_name = "longterm-sumrate"; plot_val = "throughputs"
fig = PyPlot.figure()
ax = fig[:add_axes]((0.11,0.15,0.97-0.11,0.96-0.15))

lines = Any[]
for method in [:BranchAndBoundClustering, :CoalitionFormationClustering_AttachOrSupplant, :Peters2012_Heuristic, :NoClustering, :GrandCoalitionClustering]
    line = ax[:plot](xvals, yvals_assignment(string(method), plot_val),
        color=colours_assignment[method],
        linestyle=linestyles_assignment[method],
        marker=markers_assignment[method], markeredgecolor=colours_assignment[method], markevery=15,
        label=labels_assignment[method])
    push!(lines, line[1])
end
ax[:set_ylim]([-1, 40])
ax[:set_xlabel](L"Frame structure parameter $\beta$")
ax[:set_ylabel]("IIA sum throughput [bits/s/Hz]")
legend1 = PyPlot.legend(handles=lines[1:3], loc="upper left")
legend1_frame = legend1[:get_frame]()
PyPlot.setp(legend1_frame, linewidth=0.5)
ax[:add_artist](legend1)
legend2 = ax[:legend](handles=lines[4:5], bbox_to_anchor=[0.5, 0.15], loc="center")
legend2_frame = legend2[:get_frame]()
PyPlot.setp(legend2_frame, linewidth=0.5)
fig[:savefig]("$(sim_name)_$(plot_name).eps")

##########################################################################
# longterm-sumrate_split
plot_name = "longterm-sumrate_split"
fig = PyPlot.figure()
ax = fig[:add_axes]((0.11,0.15,0.97-0.11,0.96-0.15))

for method in [:BranchAndBoundClustering, :CoalitionFormationClustering_AttachOrSupplant, :Peters2012_Heuristic]
    ax[:plot](xvals, yvals_assignment(string(method), "throughputs_network_sdma"),
        color=colours_assignment[method],
        linestyle="-",
        marker=markers_assignment[method], markeredgecolor=colours_assignment[method], markevery=15,
        label=labels_assignment[method])
end
for method in [:BranchAndBoundClustering, :CoalitionFormationClustering_AttachOrSupplant]
    ax[:plot](xvals, yvals_assignment(string(method), "throughputs_cluster_sdma"),
        color=colours_assignment[method],
        linestyle="--",
        marker=markers_assignment[method], markeredgecolor=colours_assignment[method], markevery=15)
end
ax[:set_ylim]([-1, 40])
ax[:set_xlabel](L"Frame structure parameter $\beta$")
ax[:set_ylabel]("IIA sum throughput [bits/s/Hz]")
show_legend!(ax, "upper left")
fig[:savefig]("$(sim_name)_$(plot_name).eps")
