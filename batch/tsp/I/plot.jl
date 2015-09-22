#!/usr/bin/env julia

include(joinpath(dirname(@__FILE__), "../../../src/DistributedBSClustering.jl"))
using CoordinatedPrecoding, DistributedBSClustering
using Compat, JLD, LaTeXStrings

include(joinpath(dirname(@__FILE__), "../simulation_params.jl"))
include(joinpath(dirname(@__FILE__), "../simulation_params-I.jl"))
include(joinpath(dirname(@__FILE__), "../plot_params-final.jl"))

##########################################################################
# Setup
sim_name = "I"
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

##########################################################################
# longterm-num_searches
plot_name = "longterm-num_searches"
fig = PyPlot.figure()
ax = fig[:add_axes]((0.12,0.15,0.98-0.12,0.96-0.15))

num_searches1 = mean(data["results_searches"][:,:,1], 2)
num_searches2 = mean(data["results_searches"][:,:,2], 2)
throughputs1 = mean(data["results_throughputs"][:,:,1], 2)
throughputs2 = mean(data["results_throughputs"][:,:,2], 2)

ax[:plot](data["Is"], num_searches1,
    color=colours_assignment[:CoalitionFormationClustering_AttachOrSupplant], linestyle="-", marker="o", markeredgecolor=colours_assignment[:CoalitionFormationClustering_AttachOrSupplant], markevery=2,
    label=labels_assignment[:CoalitionFormationClustering_AttachOrSupplant])
ax[:plot](data["Is"], num_searches2,
    color=colours_assignment[:CoalitionFormationClustering_Attach], linestyle="-", marker="o", markeredgecolor=colours_assignment[:CoalitionFormationClustering_Attach], markevery=2,
    label=labels_assignment[:CoalitionFormationClustering_Attach])

for (idx, It) in enumerate(data["Is"])
    if mod(idx, 2) == 1
        ax[:text](It, num_searches1[idx]-0.3, @sprintf("%.1f", throughputs1[idx]), fontdict=@compat Dict(:size => 5))
        ax[:text](It, num_searches2[idx]-0.3, @sprintf("%.1f", throughputs2[idx]), fontdict=@compat Dict(:size => 5))
    end
end

ax[:set_xlim]([0, 32])
ax[:set_ylim]([0, 8])
ax[:set_xlabel](L"Number of cells $I$")
ax[:set_ylabel](L"Average number of searches $\eta_i$")
legend = ax[:legend](loc="upper left")
legend_frame = legend[:get_frame]()
PyPlot.setp(legend_frame, linewidth=0.5)
fig[:savefig]("$(sim_name)_$(plot_name).eps")
