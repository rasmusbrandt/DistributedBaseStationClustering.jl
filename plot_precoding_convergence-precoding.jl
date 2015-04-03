#!/usr/bin/env julia

##########################################################################
# plot_precoding_convergence-precoding.jl
#
# Plots convergence curves, comparing different precoding methods.
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
        :ylim => [0, 70],
    ],

    "legend" => [
        :loc => "best",
        :fontsize => 6,
        :ncol => 4,
    ],

    "methods" => [
        "RobustIntraclusterWMMSE" => [
            ("weighted_logdet_rates", [ :color => "m", :linestyle => "-", :label => "RobustIntraclusterWMMSE" ]),
            ("utilities", [ :color => "m", :linestyle => "--", :label => "RobustIntraclusterWMMSE (utils.)" ]),
        ],

        "NaiveIntraclusterWMMSE" => [
            ("weighted_logdet_rates", [ :color => "m", :linestyle => ":", :label => "NaiveIntraclusterWMMSE" ]),
            ("utilities", [ :color => "m", :linestyle => ".", :label => "NaiveIntraclusterWMMSE (utils.)" ]),
        ],


        "RobustIntraclusterLeakageMinimization" => [
            ("weighted_logdet_rates", [ :color => "gray", :linestyle => "-", :label => "RobustIntraclusterLeakageMinimization" ]),
        ],

        "NaiveIntraclusterLeakageMinimization" => [
            ("weighted_logdet_rates", [ :color => "gray", :linestyle => ":", :label => "NaiveIntraclusterLeakageMinimization" ]),
        ],


        "RobustChen2014_MaxSINR" => [
            ("weighted_logdet_rates", [ :color => "y", :linestyle => "-", :label => "RobustChen2014_MaxSINR" ]),
        ],

        "NaiveChen2014_MaxSINR" => [
            ("weighted_logdet_rates", [ :color => "y", :linestyle => ":", :label => "NaiveChen2014_MaxSINR" ]),
        ],


        "Shi2011_WMMSE" => [
            ("weighted_logdet_rates", [ :color => "b", :linestyle => "-", :label => "WMMSE" ]),
        ],

        "Eigenprecoding" => [
            ("intercell_tdma_weighted_logdet_rates", [ :color => "c", :linestyle => "-", :label => "TDMA" ]),
            ("intracell_tdma_weighted_logdet_rates", [ :color => "c", :linestyle => "-.",  :label => "Intracell TDMA" ]),
            ("uncoord_weighted_logdet_rates", [ :color => "k", :linestyle => "-", :label => "Uncoord. transm." ]),
        ],
    ]
]

##########################################################################
# Plot it
for file_name in parsed_args["file_names"]
    data = load(file_name)
    processed_results = postprocess_precoding_convergence(data["raw_results"], data["simulation_params"], plot_params)
    plot_precoding_convergence(processed_results, data["simulation_params"], plot_params)
end
