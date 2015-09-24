#!/usr/bin/env julia

using CoordinatedPrecoding, DistributedBaseStationClustering
using Compat, JLD, LaTeXStrings

include(joinpath(dirname(@__FILE__), "../simulation_params.jl"))
include(joinpath(dirname(@__FILE__), "../simulation_params-assignment_methods.jl"))
include(joinpath(dirname(@__FILE__), "../simulation_params-num_coherence_symbols.jl"))
include(joinpath(dirname(@__FILE__), "../plot_params-final.jl"))

##########################################################################
# Setup
sim_name = "num_coherence_symbols"
data = load("$(sim_name).jld")

# Defaults
PyPlot.rc("lines", linewidth=1., markersize=3, markeredgewidth=0.5)
PyPlot.rc("font", size=7, family="serif", serif="Computer Modern Sans Serif")
PyPlot.rc("text", usetex=true)
PyPlot.rc("text.latex", preamble="\\usepackage{amsmath}")
PyPlot.rc("axes", linewidth=0.5, labelsize=7)
PyPlot.rc("xtick", labelsize=6)
PyPlot.rc("ytick", labelsize=6)
PyPlot.rc("legend", fancybox=true, fontsize=5)
PyPlot.rc("figure", figsize=(3.2,2), dpi=300)

# Legend creation helper
show_legend!(ax, loc="lower center") = begin
    legend = ax[:legend](loc=loc)
    legend_frame = legend[:get_frame]()
    PyPlot.setp(legend_frame, linewidth=0.5)
end

# Extract the values easily
xvals = vs_kmh
yvals_assignment(method, plot_val) = data["processed_assignment_results"][mean_idx][method][plot_val]
yvals_precoding(method, plot_val) = data["processed_precoding_results"][mean_idx][method][plot_val]

##########################################################################
# longterm-sumrate
plot_name = "longterm-sumrate"; plot_val = "throughputs"
fig = PyPlot.figure()
ax = fig[:add_axes]((0.13,0.15,0.98-0.13,0.96-0.15))

lines = Any[]
for method in [:BranchAndBoundClustering, :CoalitionFormationClustering_AttachOrSupplant, :CoalitionFormationClustering_Attach, :NoClustering, :GrandCoalitionClustering]
    line = ax[:plot](xvals, yvals_assignment(string(method), plot_val),
        color=colours_assignment[method],
        linestyle=linestyles_assignment[method],
        marker=markers_assignment[method], markeredgecolor=colours_assignment[method], markevery=15,
        label=labels_assignment[method])
    push!(lines, line[1])
end
ax[:set_ylim]([-1, 40])
ax[:set_xlabel]("MS speed [km/h]")
ax[:set_ylabel]("IIA sum throughput [bits/s/Hz]")
legend1 = PyPlot.legend(handles=lines[1:5], loc="upper left", bbox_to_anchor=(0.12, 0.40))
legend1_frame = legend1[:get_frame]()
PyPlot.setp(legend1_frame, linewidth=0.5)
ax[:add_artist](legend1)
# legend2 = ax[:legend](handles=lines[4:5], bbox_to_anchor=[0.5, 0.15], loc="center")
# legend2_frame = legend2[:get_frame]()
# PyPlot.setp(legend2_frame, linewidth=0.5)
fig[:savefig]("$(sim_name)_$(plot_name).eps")

##########################################################################
# longterm-avg_cluster_size
plot_name = "longterm-avg_cluster_size"; plot_val = "avg_cluster_size"
fig = PyPlot.figure()
ax = fig[:add_axes]((0.13,0.15,0.98-0.13,0.96-0.15))

for method in [:BranchAndBoundClustering, :CoalitionFormationClustering_AttachOrSupplant, :CoalitionFormationClustering_Attach, :NoClustering]
    if method == :BranchAndBoundClustering
        label = L"IIA sum throughput optimal $\mathcal{S}$"
    elseif method == :CoalitionFormationClustering_AttachOrSupplant
        label = "Coalition form. (attach-or-supplant)"
    elseif method == :CoalitionFormationClustering_Attach
        label = "Coalition form. (attach)"
    else
        label = labels_assignment[method]
    end
    ax[:plot](xvals, yvals_assignment(string(method), plot_val),
        color=colours_assignment[method],
        linestyle=linestyles_assignment[method],
        marker=markers_assignment[method], markeredgecolor=colours_assignment[method], markevery=15,
        label=label)
end
ax[:set_ylim]([0., 4.5])
ax[:set_yticks]([0, 1, 2, 3, 4])
ax[:set_xlabel]("MS speed [km/h]")
ax[:set_ylabel]("Average coalition size")
show_legend!(ax, "upper right")
fig[:savefig]("$(sim_name)_$(plot_name).eps")

##########################################################################
# instantaneous-sumrate
plot_name = "instantaneous-sumrate"; plot_val = "weighted_logdet_rates_full"
fig = PyPlot.figure()
ax = fig[:add_axes]((0.13,0.15,0.98-0.13,0.96-0.15))

for method in [:GrandCoalitionClustering, :CoalitionFormationClustering_AttachOrSupplant_IgnoreIAFeasibility, :BranchAndBoundClustering, :CoalitionFormationClustering_AttachOrSupplant, :NoClustering, ]
    if method == :BranchAndBoundClustering
        label = L"IIA sum throughput optimal $\mathcal{S}$"
    elseif method == :CoalitionFormationClustering_AttachOrSupplant
        label = "Coalition form. (attach-or-supplant)"
    elseif method == :CoalitionFormationClustering_Attach
        label = "Coalition form. (attach)"
    elseif method == :CoalitionFormationClustering_AttachOrSupplant_IgnoreIAFeasibility
        label = "Coalition form. (a-o-s, ignoring IA feas.)"
    else
        label = labels_assignment[method]
    end
    ax[:plot](xvals, yvals_precoding(string(method), plot_val),
        color=colours_assignment[method],
        linestyle=linestyles_assignment[method],
        marker=markers_assignment[method], markeredgecolor=colours_assignment[method], markevery=10,
        label=label)
end

ax[:annotate]("Loss of CSI feasibility", xy=(12, 3), xytext=(45, 20), arrowprops=(@compat Dict(:facecolor => "black", :arrowstyle => "->")), horizontalalignment="center", fontsize=6)
ax[:annotate]("", xy=(78 , 3), xytext=(55, 16), arrowprops=(@compat Dict(:facecolor => "black", :arrowstyle => "->")), horizontalalignment="left", fontsize=7)

ax[:set_ylim]([-2, 120])
ax[:set_xlabel]("MS speed [km/h]")
ax[:set_ylabel]("WMMSE sum throughp. [bits/s/Hz]")
legend = ax[:legend](loc="upper right", numpoints=1)
legend_frame = legend[:get_frame]()
PyPlot.setp(legend_frame, linewidth=0.5)
fig[:savefig]("$(sim_name)_$(plot_name).eps")
