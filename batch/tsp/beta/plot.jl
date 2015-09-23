#!/usr/bin/env julia

include(joinpath(dirname(@__FILE__), "../../../src/DistributedBaseStationClustering.jl"))
using CoordinatedPrecoding, DistributedBaseStationClustering
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
PyPlot.rc("font", size=7, family="serif", serif="Computer Modern Sans Serif")
PyPlot.rc("text", usetex=true)
PyPlot.rc("text.latex", preamble="\\usepackage{amsmath}")
PyPlot.rc("axes", linewidth=0.5, labelsize=7)
PyPlot.rc("xtick", labelsize=6)
PyPlot.rc("ytick", labelsize=6)
PyPlot.rc("legend", fancybox=true, fontsize=6)
PyPlot.rc("figure", figsize=(3.2,2), dpi=300)

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
ax = fig[:add_axes]((0.13,0.15,0.98-0.13,0.96-0.15))

lines = Any[]
for method in [:BranchAndBoundClustering, :CoalitionFormationClustering_AttachOrSupplant, :Peters2012_Heuristic]
    if method == :CoalitionFormationClustering_AttachOrSupplant
        label = "Coalition form. (attach-or-supplant)"
    else
        label = label=labels_assignment[method]
    end
    line = ax[:plot](xvals, yvals_assignment(string(method), plot_val)[:,2],
        color=colours_assignment[method],
        linestyle="-",
        marker=markers_assignment[method], markeredgecolor=colours_assignment[method], markevery=5,
        label=label)
    push!(lines, line[1])
    line = ax[:plot](xvals, yvals_assignment(string(method), plot_val)[:,1],
        color=colours_assignment[method],
        linestyle="--",
        marker=markers_assignment[method], markeredgecolor=colours_assignment[method], markevery=5,
        label=labels_assignment[method])
end
ax[:annotate](L"$\beta = 0.65$", xy=(5, 10), xytext=(0, 20), arrowprops=(@compat Dict(:facecolor => "black", :arrowstyle => "->")), horizontalalignment="center")
ax[:annotate](L"$\beta = 0$", xy=(25, 18), xytext=(30, 7), arrowprops=(@compat Dict(:facecolor => "black", :arrowstyle => "->")), horizontalalignment="center")

ax[:set_ylim]([-1, 70])
ax[:set_xlabel]("Signal-to-noise ratio [dB]")
ax[:set_ylabel]("IIA sum throughput [bits/s/Hz]")
legend1 = PyPlot.legend(handles=lines[1:3], loc="upper left")
legend1_frame = legend1[:get_frame]()
PyPlot.setp(legend1_frame, linewidth=0.5)
ax[:add_artist](legend1)
fig[:savefig]("$(sim_name)_$(plot_name).eps")

##########################################################################
# longterm-sumrate_split-bmax
plot_name = "longterm-sumrate_split-bmax"
fig = PyPlot.figure()
ax = fig[:add_axes]((0.13,0.15,0.98-0.13,0.96-0.15))

for method in [:BranchAndBoundClustering, :CoalitionFormationClustering_AttachOrSupplant, :Peters2012_Heuristic]
    ax[:plot](xvals, yvals_assignment(string(method), "throughputs_network_sdma")[:,2],
        color=colours_assignment[method],
        linestyle="-",
        marker=markers_assignment[method], markeredgecolor=colours_assignment[method], markevery=5,
        label=labels_assignment[method])
end
for method in [:BranchAndBoundClustering, :CoalitionFormationClustering_AttachOrSupplant, :Peters2012_Heuristic]
    ax[:plot](xvals, yvals_assignment(string(method), "throughputs_cluster_sdma")[:,2],
        color=colours_assignment[method],
        linestyle="--",
        marker=markers_assignment[method], markeredgecolor=colours_assignment[method], markevery=5)
end
# ax[:set_ylim]([-1, 40])
ax[:set_xlabel]("Signal-to-noise ratio [dB]")
ax[:set_ylabel]("IIA sum throughput [bits/s/Hz]")
show_legend!(ax, "upper left")
fig[:savefig]("$(sim_name)_$(plot_name).eps")

##########################################################################
# longterm-sumrate_split-bmin
plot_name = "longterm-sumrate_split-bmin"
fig = PyPlot.figure()
ax = fig[:add_axes]((0.13,0.15,0.98-0.13,0.96-0.15))

for method in [:BranchAndBoundClustering, :CoalitionFormationClustering_AttachOrSupplant, :Peters2012_Heuristic]
    ax[:plot](xvals, yvals_assignment(string(method), "throughputs_network_sdma")[:,1],
        color=colours_assignment[method],
        linestyle="-",
        marker=markers_assignment[method], markeredgecolor=colours_assignment[method], markevery=5,
        label=labels_assignment[method])
end
for method in [:BranchAndBoundClustering, :CoalitionFormationClustering_AttachOrSupplant, :Peters2012_Heuristic]
    ax[:plot](xvals, yvals_assignment(string(method), "throughputs_cluster_sdma")[:,1],
        color=colours_assignment[method],
        linestyle="--",
        marker=markers_assignment[method], markeredgecolor=colours_assignment[method], markevery=5)
end
# ax[:set_ylim]([-1, 40])
ax[:set_xlabel]("Signal-to-noise ratio [dB]")
ax[:set_ylabel]("IIA sum throughput [bits/s/Hz]")
show_legend!(ax, "upper left")
fig[:savefig]("$(sim_name)_$(plot_name).eps")
