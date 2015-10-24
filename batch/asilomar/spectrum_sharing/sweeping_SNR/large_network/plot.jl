#!/usr/bin/env julia

include(joinpath(dirname(@__FILE__), "../../../../../src/DistributedBaseStationClustering.jl"))
using DistributedBaseStationClustering, CoordinatedPrecoding
using Compat, JLD
using LaTeXStrings

include(joinpath(dirname(@__FILE__), "../../../simulation_params.jl"))
include(joinpath(dirname(@__FILE__), "../../../simulation_params-assignment_methods.jl"))
include(joinpath(dirname(@__FILE__), "../../../simulation_params-large_network.jl"))
include(joinpath(dirname(@__FILE__), "../../../simulation_params-SNR.jl"))
include(joinpath(dirname(@__FILE__), "../../../plot_params-assignment_methods.jl"))
include(joinpath(dirname(@__FILE__), "../../../plot_params-final.jl"))

##########################################################################
# Load data
data = load("merged.jld")
simulation_params = data["simulation_params"]
results_assignment = data["results_assignment"]
results_assignment_mean = data["results_assignment_mean"]
results_assignment_var = data["results_assignment_var"]
results_precoding = data["results_precoding"]
results_precoding_mean = data["results_precoding_mean"]
results_precoding_var = data["results_precoding_var"]
data = []

##########################################################################
# Figure properties
PyPlot.rc("lines", linewidth=1, markersize=3, markeredgewidth=0.5)
PyPlot.rc("font", size=8, family="serif", serif="Computer Modern Sans Serif")
PyPlot.rc("text", usetex=true)
PyPlot.rc("text.latex", preamble="\\usepackage{amsmath}")
PyPlot.rc("axes", linewidth=0.5, labelsize=7)
PyPlot.rc("xtick", labelsize=7)
PyPlot.rc("ytick", labelsize=7)
PyPlot.rc("legend", fancybox=true, fontsize=6, numpoints=1)
PyPlot.rc("figure", figsize=(3.5,2.0))

##########################################################################
# Plots
fig = PyPlot.figure()
ax = fig[:add_axes]((0.11,0.16,0.95-0.11,0.95-0.16))

ax[:plot](transmit_powers_dBm, results_precoding_mean["CoalitionFormationClustering_Individual"]["weighted_logdet_rates_LB"][:,2], color=colours[:CoalitionFormationClustering_Individual], markeredgecolor=colours[:CoalitionFormationClustering_Individual], linestyle="-", marker=markers[:CoalitionFormationClustering_Individual], label=L"Coalition formation ($b_k = 10$)")
ax[:plot](transmit_powers_dBm, results_precoding_mean["CoalitionFormationClustering_Individual"]["weighted_logdet_rates_LB"][:,1], color=colours[:CoalitionFormationClustering_Individual], markeredgecolor=colours[:CoalitionFormationClustering_Individual], linestyle="--", marker=markers[:CoalitionFormationClustering_Individual], label=L"Coalition formation ($b_k = 2$)")
ax[:plot](transmit_powers_dBm, results_precoding_mean["Chen2014_kmeans"]["weighted_logdet_rates_LB"][:,1], markeredgewidth=1, color=colours[:Chen2014_kmeans], markeredgecolor=colours[:Chen2014_kmeans], linestyle="-", marker=markers[:Chen2014_kmeans], label=labels[:Chen2014_kmeans])
ax[:plot](transmit_powers_dBm, results_precoding_mean["GrandCoalitionClustering"]["weighted_logdet_rates_LB"][:,1], color=colours[:GrandCoalitionClustering], markeredgecolor=colours[:GrandCoalitionClustering], linestyle="-", marker=markers[:GrandCoalitionClustering], label=labels[:GrandCoalitionClustering])
ax[:plot](transmit_powers_dBm, results_precoding_mean["NoClustering"]["weighted_logdet_rates_LB"][:,1], color=colours[:NoClustering], markeredgecolor=colours[:NoClustering], linestyle="-", marker=markers[:NoClustering], label=labels[:NoClustering])

ax[:set_ylim]([0, 65])

ax[:set_xlabel]("Transmit power [dBm]")
ax[:set_ylabel]("Instant. sum throughput [bits/s/Hz]")

legend = ax[:legend](loc="lower right")
# legend_lines = legend[:get_lines]()
legend_frame = legend[:get_frame]()
# PyPlot.setp(legend_lines, linewidth=0.5)
PyPlot.setp(legend_frame, linewidth=0.5)

fig[:savefig]("large_network-SNR-instantaneous_sumrate.eps")

PyPlot.rc("figure", figsize=(3.5,2.2))
fig = PyPlot.figure()
ax = fig[:add_subplot](2, 1, 1)

ax[:plot](transmit_powers_dBm, simulation_params["I"]./results_assignment_mean["Chen2014_kmeans"]["num_clusters"][:,1], markeredgewidth=1, color=colours[:Chen2014_kmeans], markeredgecolor=colours[:Chen2014_kmeans], linestyle="-", marker=markers[:Chen2014_kmeans], label=labels[:Chen2014_kmeans])
ax[:plot](transmit_powers_dBm, simulation_params["I"]./results_assignment_mean["CoalitionFormationClustering_Individual"]["num_clusters"][:,2], color=colours[:CoalitionFormationClustering_Individual], markeredgecolor=colours[:CoalitionFormationClustering_Individual], linestyle="-", marker=markers[:CoalitionFormationClustering_Individual], label=L"Coalition formation ($b_k = 10$)")
ax[:plot](transmit_powers_dBm, simulation_params["I"]./results_assignment_mean["CoalitionFormationClustering_Individual"]["num_clusters"][:,1], color=colours[:CoalitionFormationClustering_Individual], markeredgecolor=colours[:CoalitionFormationClustering_Individual], linestyle="--", marker=markers[:CoalitionFormationClustering_Individual], label=L"Coalition formation ($b_k = 2$)")

ax[:set_ylim]([0, 4.5])
ax[:set_yticks](0:4)

ax[:set_ylabel](L"Coalition size $\lvert \mathcal{C}_p \rvert$")

legend = ax[:legend](loc="lower right")
# legend_lines = legend[:get_lines]()
legend_frame = legend[:get_frame]()
# PyPlot.setp(legend_lines, linewidth=0.5)
PyPlot.setp(legend_frame, linewidth=0.5)


ax = fig[:add_subplot](2, 1, 2)

ax[:plot](transmit_powers_dBm, (1/simulation_params["I"])*results_assignment_mean["CoalitionFormationClustering_Individual"]["num_searches"][:,2], color=colours[:CoalitionFormationClustering_Individual], markeredgecolor=colours[:CoalitionFormationClustering_Individual], linestyle="-", marker=markers[:CoalitionFormationClustering_Individual])
ax[:plot](transmit_powers_dBm, (1/simulation_params["I"])*results_assignment_mean["CoalitionFormationClustering_Individual"]["num_searches"][:,1], color=colours[:CoalitionFormationClustering_Individual], markeredgecolor=colours[:CoalitionFormationClustering_Individual], linestyle="--", marker=markers[:CoalitionFormationClustering_Individual])

ax[:set_ylim]([0, 4])
ax[:set_yticks](0:4)

ax[:set_xlabel]("Transmit power [dBm]")
ax[:set_ylabel](L"\# of deviations $\eta_k$")

fig[:tight_layout]()

fig[:savefig]("large_network-SNR-num_clusters-num_searches.eps")
