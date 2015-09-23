#!/usr/bin/env julia

include(joinpath(dirname(@__FILE__), "../../../src/DistributedBaseStationClustering.jl"))
using CoordinatedPrecoding, DistributedBaseStationClustering
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
PyPlot.rc("font", size=7, family="serif", serif="Computer Modern Sans Serif")
PyPlot.rc("text", usetex=true)
PyPlot.rc("text.latex", preamble="\\usepackage{amsmath}")
PyPlot.rc("axes", linewidth=0.5, labelsize=7)
PyPlot.rc("xtick", labelsize=6)
PyPlot.rc("ytick", labelsize=6)
PyPlot.rc("legend", fancybox=true, fontsize=6)
PyPlot.rc("figure", figsize=(3.2,2), dpi=300)

##########################################################################
# longterm-num_searches
plot_name = "longterm-num_searches"
fig = PyPlot.figure()
ax = fig[:add_axes]((0.13,0.15,0.98-0.13,0.96-0.15))

num_searches1 = mean(data["results_searches"][:,:,1], 2)
num_searches2 = mean(data["results_searches"][:,:,2], 2)
throughputs1 = mean(data["results_throughputs"][:,:,1], 2)
throughputs2 = mean(data["results_throughputs"][:,:,2], 2)

ax[:plot](data["Is"], num_searches1,
    color=colours_assignment[:CoalitionFormationClustering_AttachOrSupplant], linestyle="-", marker="o", markeredgecolor=colours_assignment[:CoalitionFormationClustering_AttachOrSupplant], markevery=1,
    label=labels_assignment[:CoalitionFormationClustering_AttachOrSupplant])
ax[:plot](data["Is"], num_searches2,
    color=colours_assignment[:CoalitionFormationClustering_Attach], linestyle="-", marker="o", markeredgecolor=colours_assignment[:CoalitionFormationClustering_Attach], markevery=1,
    label=labels_assignment[:CoalitionFormationClustering_Attach])
ax[:plot](1:9, 0:8, "k:") # -1 since no searches are needed for I = 1
ax[:annotate]("Linear growth", xy=(7, 5), xytext=(17, 7), arrowprops=(@compat Dict(:facecolor => "black", :arrowstyle => "->")), horizontalalignment="center", fontsize=6)

for (idx, It) in enumerate(data["Is"][1:end-1])
    if idx == 1; continue; end
    ax[:text](It, num_searches1[idx]-0.5, @sprintf("%.1f", throughputs1[idx]), fontdict=@compat Dict(:size => 5))
    ax[:text](It, num_searches2[idx]-0.5, @sprintf("%.1f", throughputs2[idx]), fontdict=@compat Dict(:size => 5))
end

ax[:set_xlim]([0, 50])
ax[:set_xticks]([2, 10, 20, 30, 40, 50])
ax[:set_ylim]([0, 8])
ax[:set_xlabel](L"Number of cells $I$")
ax[:set_ylabel](L"Average number of searches $\eta_i$")
legend = ax[:legend](loc="upper left", bbox_to_anchor=(0.18, 0.4))
legend_frame = legend[:get_frame]()
PyPlot.setp(legend_frame, linewidth=0.5)
fig[:savefig]("$(sim_name)_$(plot_name).eps")
