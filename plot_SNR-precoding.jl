#!/usr/bin/env julia

##########################################################################
# plot_system-SNR-precoding.jl
#
# Plots SNR curves, comparing different precoding methods.
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
        :fontsize => 6,
    ],

    "methods" => [
        "RobustIntraclusterWMMSE" => [
            ("weighted_logdet_rates_full", [ :color => "m", :linestyle => "-", :label => "RobustIntraclusterWMMSE (full CSI)" ]),
            ("weighted_logdet_rates_partial", [ :color => "m", :linestyle => "--", :label => "RobustIntraclusterWMMSE (partial CSI)" ]),
        ],

        "NaiveIntraclusterWMMSE" => [
            ("weighted_logdet_rates_full", [ :color => "c", :linestyle => "-", :label => "NaiveIntraclusterWMMSE (full CSI)" ]),
            ("weighted_logdet_rates_partial", [ :color => "c", :linestyle => "--", :label => "NaiveIntraclusterWMMSE (partial CSI)" ]),
        ],


        "RobustIntraclusterLeakageMinimization" => [
            ("weighted_logdet_rates_full", [ :color => "DarkGray", :linestyle => "-", :label => "RobustIntraclusterLeakageMinimization (full CSI)" ]),
            ("weighted_logdet_rates_partial", [ :color => "DarkGray", :linestyle => "--", :label => "RobustIntraclusterLeakageMinimization (partial CSI)" ]),
        ],

        "NaiveIntraclusterLeakageMinimization" => [
            ("weighted_logdet_rates_full", [ :color => "gray", :linestyle => "-", :label => "NaiveIntraclusterLeakageMinimization (full CSI)" ]),
            ("weighted_logdet_rates_partial", [ :color => "gray", :linestyle => "--", :label => "NaiveIntraclusterLeakageMinimization (partial CSI)" ]),
        ],


        "RobustChen2014_MaxSINR" => [
            ("weighted_logdet_rates_full", [ :color => "g", :linestyle => "-", :label => "RobustChen2014_MaxSINR (full CSI)" ]),
            ("weighted_logdet_rates_partial", [ :color => "g", :linestyle => "--", :label => "RobustChen2014_MaxSINR (partial CSI)" ]),
        ],

        "NaiveChen2014_MaxSINR" => [
            ("weighted_logdet_rates_full", [ :color => "y", :linestyle => "-", :label => "NaiveChen2014_MaxSINR (full CSI)" ]),
            ("weighted_logdet_rates_partial", [ :color => "y", :linestyle => "--", :label => "NaiveChen2014_MaxSINR (partial CSI)" ]),
        ],


        "Shi2011_WMMSE" => [
            ("weighted_logdet_rates", [ :color => "b", :linestyle => "-", :label => "WMMSE (full CSI)" ]),
        ],

        "Eigenprecoding" => [
            ("intercell_tdma_weighted_logdet_rates", [ :color => "k", :linestyle => "-", :label => "TDMA" ]),
            ("intracell_tdma_weighted_logdet_rates", [ :color => "k", :linestyle => "--",  :label => "Intracell TDMA" ]),
            ("uncoord_weighted_logdet_rates", [ :color => "k", :linestyle => ":", :label => "Uncoord. transm." ]),
        ],
    ]
]

##########################################################################
# Plot it
for file_name in parsed_args["file_names"]
    data = load(file_name)
    processed_results = postprocess(data["raw_precoding_results"], data["simulation_params"], plot_params)
    plot(processed_results, data["simulation_params"], plot_params)
end
