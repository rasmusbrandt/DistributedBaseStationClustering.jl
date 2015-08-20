#!/usr/bin/env julia

include(joinpath(dirname(@__FILE__), "../../../src/IAClustering.jl"))
using CoordinatedPrecoding, IAClustering
using JLD, LaTeXStrings

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
PyPlot.rc("legend", fancybox=true, fontsize=8)
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
ax = fig[:add_axes]((0.11,0.15,0.95-0.11,0.95-0.15))

for method in [:BranchAndBoundClustering, :CoalitionFormationClustering_AttachOrSupplant, :NoClustering, :GrandCoalitionClustering]
    if method == :CoalitionFormationClustering_AttachOrSupplant
        ax[:plot](xvals, yvals_assignment(string(method), plot_val)[:,2],
            color=colours_assignment[method],
            linestyle="-",
            marker=markers_assignment[method], markeredgecolor=colours_assignment[method], markevery=15,
            label=latexstring("$(labels_assignment[method]) (\$b_i = \\infty\$)"))
        ax[:plot](xvals, yvals_assignment(string(method), plot_val)[:,1],
            color=colours_assignment[method],
            linestyle="--",
            marker=markers_assignment[method], markeredgecolor=colours_assignment[method], markevery=15,
            label=latexstring("$(labels_assignment[method]) (\$b_i = 5\$)"))
    else
        ax[:plot](xvals, yvals_assignment(string(method), plot_val)[:,1],
            color=colours_assignment[method],
            linestyle=linestyles_assignment[method],
            marker=markers_assignment[method], markeredgecolor=colours_assignment[method], markevery=15,
            label=labels_assignment[method])
    end
end
ax[:set_ylim]([-1, 60])
ax[:set_xlabel]("Signal-to-noise ratio [dB]")
ax[:set_ylabel]("Long-term sum throughput [bits/s/Hz]")
show_legend!(ax, "upper left")
fig[:savefig]("$(sim_name)_$(plot_name).eps")

##########################################################################
# longterm-num_searches
plot_name = "longterm-num_searches"; plot_val = "num_searches"
fig = PyPlot.figure()
ax = fig[:add_axes]((0.11,0.15,0.95-0.11,0.95-0.15))

for method in [:CoalitionFormationClustering_AttachOrSupplant]
    if method == :CoalitionFormationClustering_AttachOrSupplant
        ax[:plot](xvals, (1/I)*yvals_assignment(string(method), plot_val)[:,2],
            color=colours_assignment[method],
            linestyle="-",
            marker=markers_assignment[method], markeredgecolor=colours_assignment[method], markevery=15,
            label=latexstring("$(labels_assignment[method]) (\$b_i = \\infty\$)"))
        ax[:plot](xvals, (1/I)*yvals_assignment(string(method), plot_val)[:,1],
            color=colours_assignment[method],
            linestyle="--",
            marker=markers_assignment[method], markeredgecolor=colours_assignment[method], markevery=15,
            label=latexstring("$(labels_assignment[method]) (\$b_i = 5\$)"))
    else
        ax[:plot](xvals, (1/I)*yvals_assignment(string(method), plot_val)[:,1],
            color=colours_assignment[method],
            linestyle=linestyles_assignment[method],
            marker=markers_assignment[method], markeredgecolor=colours_assignment[method], markevery=15,
            label=labels_assignment[method])
    end
end
ax[:set_xlabel]("Signal-to-noise ratio [dB]")
ax[:set_ylabel](L"Average number of searches $\eta_i$")
show_legend!(ax, "upper left")
fig[:savefig]("$(sim_name)_$(plot_name).eps")

##########################################################################
# instantaneous-sumrate
plot_name = "instantaneous-sumrate"; plot_val = "weighted_logdet_rates_full"
fig = PyPlot.figure()
ax = fig[:add_axes]((0.11,0.15,0.95-0.11,0.95-0.15))

for method in [:BranchAndBoundClustering, :CoalitionFormationClustering_AttachOrSupplant, :Chen2014_kmeans, :NoClustering, :GrandCoalitionClustering]
    if method == :CoalitionFormationClustering_AttachOrSupplant
        ax[:plot](xvals, yvals_precoding(string(method), plot_val)[:,2],
            color=colours_assignment[method],
            linestyle="-",
            marker=markers_assignment[method], markeredgecolor=colours_assignment[method], markevery=15,
            label=latexstring("$(labels_assignment[method]) (\$b_i = \\infty\$)"))
        ax[:plot](xvals, yvals_precoding(string(method), plot_val)[:,1],
            color=colours_assignment[method],
            linestyle="--",
            marker=markers_assignment[method], markeredgecolor=colours_assignment[method], markevery=15,
            label=latexstring("$(labels_assignment[method]) (\$b_i = 5\$)"))
    else
        ax[:plot](xvals, yvals_precoding(string(method), plot_val)[:,1],
            color=colours_assignment[method],
            linestyle=linestyles_assignment[method],
            marker=markers_assignment[method], markeredgecolor=colours_assignment[method], markevery=15,
            label=labels_assignment[method])
    end
end
ax[:set_ylim]([-2, 100])
ax[:set_xlabel]("Signal-to-noise ratio [dB]")
ax[:set_ylabel]("Average sum throughput [bits/s/Hz]")
show_legend!(ax, "upper left")
fig[:savefig]("$(sim_name)_$(plot_name).eps")
