#!/usr/bin/env julia

include(joinpath(dirname(@__FILE__), "../../../../../src/DistributedBaseStationClustering.jl"))
using DistributedBaseStationClustering, CoordinatedPrecoding
using Compat, JLD
using LaTeXStrings

include(joinpath(dirname(@__FILE__), "../../../simulation_params.jl"))
include(joinpath(dirname(@__FILE__), "../../../simulation_params-assignment_methods.jl"))
include(joinpath(dirname(@__FILE__), "../../../simulation_params-small_network1.jl"))
include(joinpath(dirname(@__FILE__), "../../../simulation_params-num_coherence_symbols.jl"))
include(joinpath(dirname(@__FILE__), "../../../plot_params-assignment_methods.jl"))
include(joinpath(dirname(@__FILE__), "../../../plot_params-final.jl"))

##########################################################################
# Load data
using Compat, JLD
data = load("raw-small_network1.jld")

##########################################################################
# Perform post processing
results_assignment, results_assignment_mean, results_assignment_var = postprocess(data["raw_assignment_results"], data["simulation_params"], postprocess_params_assignment)
results_precoding, results_precoding_mean, results_precoding_var = postprocess(data["raw_precoding_results"], data["simulation_params"], postprocess_params_precoding)

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

ax[:plot](vs_kmh, results_assignment_mean["ExhaustiveSearchClustering"]["utilities"], color=colours[:ExhaustiveSearchClustering], markeredgecolor=colours[:ExhaustiveSearchClustering], linestyle="-", marker=markers[:ExhaustiveSearchClustering], label=labels[:ExhaustiveSearchClustering])
ax[:plot](vs_kmh, results_assignment_mean["CoalitionFormationClustering_Individual"]["utilities"], color=colours[:CoalitionFormationClustering_Individual], markeredgecolor=colours[:CoalitionFormationClustering_Individual], linestyle="-", marker=markers[:CoalitionFormationClustering_Individual], label=L"Coalition formation ($b_k = 10$)")
ax[:plot](vs_kmh, results_assignment_mean["GrandCoalitionClustering"]["utilities"], color=colours[:GrandCoalitionClustering], markeredgecolor=colours[:GrandCoalitionClustering], linestyle="-", marker=markers[:GrandCoalitionClustering], label=labels[:GrandCoalitionClustering])
ax[:plot](vs_kmh, results_assignment_mean["RandomClustering"]["utilities"], color=colours[:RandomClustering], markeredgecolor=colours[:RandomClustering], linestyle="-", marker=markers[:RandomClustering], label=labels[:RandomClustering])
ax[:plot](vs_kmh, results_assignment_mean["NoClustering"]["utilities"], color=colours[:NoClustering], markeredgecolor=colours[:NoClustering], linestyle="-", marker=markers[:NoClustering], label=labels[:NoClustering])

ax[:set_ylim]([-1, 26])

ax[:set_xlabel]("MS speed [km/h]")
ax[:set_ylabel]("Long-term sum throughput [bits/s/Hz]", fontsize=6)

legend = ax[:legend](loc="best")
# legend_lines = legend[:get_lines]()
legend_frame = legend[:get_frame]()
# PyPlot.setp(legend_lines, linewidth=0.5)
PyPlot.setp(legend_frame, linewidth=0.5)

fig[:savefig]("small_network1-coherence-longterm_sumrate.eps")
