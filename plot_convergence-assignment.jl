#!/usr/bin/env julia

##########################################################################
# plot_convergence-cluster_assignment.jl
#
# Plots convergence curves, comparing different cluster assignment methods.
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

    "objective" => :sumrate,

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
    ],

    "methods" => [
        "GrandCoalitionClustering" => [
            ("logdet_rates", [ :color => "b", :linestyle => "-", :label => "GrandCoalitionClustering" ]),
        ],

        "ExhaustiveSearchClustering" => [
            ("logdet_rates", [ :color => "y", :linestyle => "-", :label => "ExhaustiveSearchClustering" ]),
            ("utilities", [ :color => "y", :linestyle => "--", :label => "ExhaustiveSearchClustering (utils.)" ]),
        ],

        "BnBClustering" => [
            ("logdet_rates", [ :color => "m", :linestyle => "-", :label => "BnBClustering" ]),
            ("utilities", [ :color => "m", :linestyle => "--", :label => "BnBClustering (utils.)" ]),
        ],

        "Chen2014_LinearObjClustering_Exhaustive" => [
            ("logdet_rates", [ :color => "g", :linestyle => "-", :label => "Chen2014_LinearObjClustering_Exhaustive" ]),
            ("utilities", [ :color => "g", :linestyle => "--", :label => "Chen2014_LinearObjClustering_Exhaustive (utils.)" ]),
        ],

        "NeighbourClustering" => [
            ("logdet_rates", [ :color => "c", :linestyle => "-", :label => "NeighbourClustering" ]),
            ("utilities", [ :color => "c", :linestyle => "--", :label => "NeighbourClustering (utils.)" ]),
        ],

        "RandomClustering" => [
            ("logdet_rates", [ :color => "k", :linestyle => "-", :label => "RandomClustering" ]),
            ("utilities", [ :color => "k", :linestyle => "--", :label => "RandomClustering (utils.)" ]),
        ],
    ]
]

##########################################################################
# Plot it
for file_name in parsed_args["file_names"]
    data = load(file_name)
    processed_results = postprocess_convergence(data["raw_results"], data["simulation_params"], plot_params)
    plot_convergence(processed_results, data["simulation_params"], plot_params)
end
