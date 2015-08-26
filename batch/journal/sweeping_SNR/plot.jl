#!/usr/bin/env julia

include(joinpath(dirname(@__FILE__), "../../../src/IAClustering.jl"))
using CoordinatedPrecoding, IAClustering
using Compat, JLD, LaTeXStrings

include(joinpath(dirname(@__FILE__), "../simulation_params.jl"))
include(joinpath(dirname(@__FILE__), "../simulation_params-assignment_methods.jl"))
include(joinpath(dirname(@__FILE__), "../simulation_params-SNR.jl"))
include(joinpath(dirname(@__FILE__), "../plot_params-final.jl"))

##########################################################################
# Setup
sim_name = "SNR"
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
ax = fig[:add_axes]((0.11,0.15,0.98-0.11,0.96-0.15))

lines = Any[]
for method in [:BranchAndBoundClustering, :CoalitionFormationClustering_AttachOrSupplant, :Peters2012_Heuristic, :NoClustering, :GrandCoalitionClustering]
    if method == :CoalitionFormationClustering_AttachOrSupplant
        line = ax[:plot](xvals, yvals_assignment(string(method), plot_val)[:,2],
            color=colours_assignment[method],
            linestyle="-",
            marker=markers_assignment[method], markeredgecolor=colours_assignment[method], markevery=5,
            label=latexstring("Coalition form. (\$b_i = \\infty\$)"))
        line = ax[:plot](xvals, yvals_assignment(string(method), plot_val)[:,1],
            color=colours_assignment[method],
            linestyle="--",
            marker=markers_assignment[method], markeredgecolor=colours_assignment[method], markevery=5,
            label=latexstring("Coalition form. (\$b_i = 5\$)"))
    else
        line = ax[:plot](xvals, yvals_assignment(string(method), plot_val)[:,1],
            color=colours_assignment[method],
            linestyle=linestyles_assignment[method],
            marker=markers_assignment[method], markeredgecolor=colours_assignment[method], markevery=5,
            label=labels_assignment[method])
    end
    push!(lines, line[1])
end
ax[:set_ylim]([-1, 60])
ax[:set_xlabel]("Signal-to-noise ratio [dB]")
ax[:set_ylabel]("IIA sum throughput [bits/s/Hz]")
show_legend!(ax, "upper left")
fig[:savefig]("$(sim_name)_$(plot_name).eps")

##########################################################################
# longterm-num_searches
plot_name = "longterm-num_searches"; plot_val = "num_searches"
fig = PyPlot.figure()
ax = fig[:add_axes]((0.09,0.15,0.98-0.09,0.96-0.15))

method = :CoalitionFormationClustering_AttachOrSupplant
ax[:plot](xvals, (1/I)*yvals_assignment(string(method), plot_val)[:,2],
    color=colours_assignment[method],
    linestyle="-",
    marker=markers_assignment[method], markeredgecolor=colours_assignment[method], markevery=5,
    label=latexstring("$(labels_assignment[method]) (\$b_i = \\infty\$)"))
ax[:plot](xvals, (1/I)*yvals_assignment(string(method), plot_val)[:,1],
    color=colours_assignment[method],
    linestyle="--",
    marker=markers_assignment[method], markeredgecolor=colours_assignment[method], markevery=5,
    label=latexstring("$(labels_assignment[method]) (\$b_i = 5\$)"))
ax[:set_xlabel]("Signal-to-noise ratio [dB]")
ax[:set_ylabel](L"Average number of searches $\eta_i$")
show_legend!(ax, "upper right")
fig[:savefig]("$(sim_name)_$(plot_name).eps")

##########################################################################
# instantaneous-sumrate
plot_name = "instantaneous-sumrate"; plot_val = "weighted_logdet_rates_full"
fig = PyPlot.figure()
ax = fig[:add_axes]((0.13,0.15,0.98-0.13,0.96-0.15))

lines = Any[]
for method in [:BranchAndBoundClustering, :CoalitionFormationClustering_AttachOrSupplant, :Peters2012_Heuristic, :NoClustering, :Chen2014_kmeans, :GrandCoalitionClustering]
    if method == :CoalitionFormationClustering_AttachOrSupplant
        line1 = ax[:plot](xvals, yvals_precoding(string(method), plot_val)[:,2],
            color=colours_assignment[method],
            linestyle="-",
            marker=markers_assignment[method], markeredgecolor=colours_assignment[method], markevery=5,
            label=latexstring("$(labels_assignment[method]) (\$b_i = \\infty\$)"))
        line2 = ax[:plot](xvals, yvals_precoding(string(method), plot_val)[:,1],
            color=colours_assignment[method],
            linestyle="--",
            marker=markers_assignment[method], markeredgecolor=colours_assignment[method], markevery=5,
            label=latexstring("$(labels_assignment[method]) (\$b_i = 5\$)"))
        push!(lines, line1[1])
        push!(lines, line2[1])
    else
        line = ax[:plot](xvals, yvals_precoding(string(method), plot_val)[:,1],
            color=colours_assignment[method],
            linestyle=linestyles_assignment[method],
            marker=markers_assignment[method], markeredgecolor=colours_assignment[method], markevery=5,
            label=labels_assignment[method])
        push!(lines, line[1])
    end
end
ax[:set_ylim]([-2, 100])
ax[:set_xlabel]("Signal-to-noise ratio [dB]")
ax[:set_ylabel]("WMMSE sum throughput [bits/s/Hz]")
legend1 = PyPlot.legend(handles=lines[1:3], loc="upper left")
legend1_frame = legend1[:get_frame]()
PyPlot.setp(legend1_frame, linewidth=0.5)
ax[:add_artist](legend1)
legend2 = ax[:legend](handles=lines[4:7], bbox_to_anchor=[0.5, 0.05], loc="lower left")
legend2_frame = legend2[:get_frame]()
PyPlot.setp(legend2_frame, linewidth=0.5)
fig[:savefig]("$(sim_name)_$(plot_name).eps")
