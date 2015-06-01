#!/usr/bin/env julia

##########################################################################
# plot_assignment_convergence.jl
#
# Plots convergence curves for assignment methods.
##########################################################################

include("src/IAClustering.jl")
using IAClustering, CoordinatedPrecoding

##########################################################################
# Load data
#
# Do this before loading other code, otherwise the JLD module might crash!
using HDF5, JLD, ArgParse
s = ArgParseSettings()
@add_arg_table s begin
    "file_names"
        help = "file names with results"
        required = true
        nargs = '+'
end
parsed_args = parse_args(s)

##########################################################################
# Plot parameters
plot_params = [
    "plot_name" => "",

    "objective" => :none,

    "figure" => [
        :figsize => (8,5),
        :dpi => 125,
    ],

    "axes" => [
        :xlabel => "Iterations",
        :ylabel => "Sum rate [bits/s/Hz]",
    ],

    "legend" => [
        :loc => "best",
        :fontsize => 8,
        :ncol => 2,
    ],

    "methods" => [
        "BranchAndBoundClustering" => [
            ("upper_bound_evolution", [ :color => "Coral", :linestyle => "--", :label => "BranchAndBoundClustering (UB)" ]),
            ("lower_bound_evolution", [ :color => "Coral", :linestyle => "-", :label => "BranchAndBoundClustering (LB)" ]),
        ],
    ]
]

##########################################################################
# Plot it
for file_name in parsed_args["file_names"]
    data = load(file_name)
    processed_results = postprocess_assignment_convergence(data["raw_results"], data["simulation_params"], plot_params)
    plot_assignment_convergence(processed_results, data["simulation_params"], plot_params)
end
