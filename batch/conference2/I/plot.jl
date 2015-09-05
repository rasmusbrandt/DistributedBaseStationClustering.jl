#!/usr/bin/env julia

include("../../../src/IAClustering.jl")
using IAClustering, CoordinatedPrecoding
using Compat, JLD, LaTeXStrings

sim_name = "I"
data = load("$(sim_name).jld")

# 8-class Set1
colours = [
    :red => "#e41a1c",
    :blue => "#377eb8",
    :green => "#4daf4a",
    :purple => "#984ea3",
    :orange => "#ff7f00",
    :yellow => "#ffff33",
    :brown => "#a65628",
    :pink => "#f781bf",
]

# Plot defaults
PyPlot.rc("lines", linewidth=1., markersize=3, markeredgewidth=0.5)
PyPlot.rc("font", size=6, family="serif", serif="Computer Modern Sans Serif")
PyPlot.rc("text", usetex=true)
PyPlot.rc("text.latex", preamble="\\usepackage{amsmath}")
PyPlot.rc("axes", linewidth=0.5, labelsize=6)
PyPlot.rc("xtick", labelsize=6)
PyPlot.rc("ytick", labelsize=6)
PyPlot.rc("legend", fancybox=true, fontsize=6)
PyPlot.rc("figure", figsize=(3.50,1.5), dpi=125)

# Plot it
fig = PyPlot.figure()
ax = fig[:add_axes]((0.12,0.18,0.90-0.12,0.95-0.18))

ax[:plot](data["Is"], mean(data["results"], 2),
    color=colours[:blue], linestyle="-",
    label="Branch and bound")
ax[:set_xlabel](L"Number of BSs $I$")
ax[:set_ylabel]("Number of nodes bounded")
ax[:set_yscale]("log")
legend = ax[:legend](loc="upper left")
legend_frame = legend[:get_frame]()
PyPlot.setp(legend_frame, linewidth=0.5)
fig[:savefig]("I.eps")
