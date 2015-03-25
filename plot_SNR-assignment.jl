#!/usr/bin/env julia

##########################################################################
# plot_SNR-cluster_assignment.jl
#
# Plots SNR curves, comparing different cluster assignment methods.
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
        :xlabel => "Transmit power [dBm]",
        :ylabel => "Sum rate [bits/s/Hz]",
        :ylim => [0, 70],
    ],

    "legend" => [
        :loc => "best",
        :fontsize => 10,
    ],

    "methods" => [
        "Chen2014_ExhaustiveSearch" => [
            ("logdet_rates", [ :color => "DodgerBlue", :linestyle => "-", :label => "Chen2014_ExhaustiveSearch" ]),
        ],

        "ExhaustiveSearchClustering" => [
            ("logdet_rates", [ :color => "Coral", :linestyle => "", :marker => ".", :label => "ExhaustiveSearchClustering" ]),
        ],

        "BranchAndBoundClustering" => [
            ("logdet_rates", [ :color => "Coral", :linestyle => "-", :label => "BranchAndBoundClustering" ]),
        ],

        "GrandCoalitionClustering" => [
            ("logdet_rates", [ :color => "Maroon", :linestyle => "-", :label => "GrandCoalitionClustering" ]),
        ],

        "GreedyClustering" => [
            ("logdet_rates", [ :color => "DarkOrchid", :linestyle => "-", :label => "GreedyClustering" ]),
        ],

        "RandomClustering" => [
            ("logdet_rates", [ :color => "Khaki", :linestyle => "-", :label => "RandomClustering" ]),
        ],

        "NoClustering" => [
            ("logdet_rates", [ :color => "Pink", :linestyle => "-", :label => "NoClustering" ]),
        ],
    ]
]

##########################################################################
# Plot it
for file_name in parsed_args["file_names"]
    data = load(file_name)
    processed_results = postprocess(data["raw_results"], data["simulation_params"], plot_params)
    plot(processed_results, data["simulation_params"], plot_params)
end
