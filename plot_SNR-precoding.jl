#!/usr/bin/env julia

##########################################################################
# plot_system-SNR-precoding.jl
#
# Plots SNR curves, comparing different precoding methods.
##########################################################################

include("src/DistributedBSClustering.jl")
using DistributedBSClustering, CoordinatedPrecoding

##########################################################################
# Load data
#
# Do this before loading other code, otherwise the JLD module might crash!
using Compat, JLD, ArgParse
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
plot_params = @compat Dict(
    "plot_name" => "",

    "objective" => :sum,

    "figure" => Dict(
        :figsize => (8,5),
        :dpi => 125,
    ),

    "axes" => Dict(
        :xlabel => "SNR [dB]",
        :ylabel => "Sum rate [bits/s/Hz]",
        :ylim => [0, 70],
    ),

    "legend" => Dict(
        :loc => "best",
        :fontsize => 6,
    ),

    "methods" => Dict(
        "RobustIntraclusterWMMSE" => [
            ("weighted_logdet_rates_full", Dict(color => "m", :linestyle => "-", :label => "RobustIntraclusterWMMSE (full CSI, sum)")),
            ("weighted_logdet_rates_partial", Dict(color => "m", :linestyle => "-", :label => "RobustIntraclusterWMMSE (partial CSI, sum)")),
            ("weighted_logdet_rates_cluster_sdma_full", Dict(color => "m", :linestyle => ":", :label => "RobustIntraclusterWMMSE (full CSI, cluster)")),
            ("weighted_logdet_rates_network_sdma_full", Dict(color => "m", :linestyle => "--", :label => "RobustIntraclusterWMMSE (full CSI, network)")),
            ("weighted_logdet_rates_network_sdma_partial", Dict(color => "m", :linestyle => "--", :label => "RobustIntraclusterWMMSE (partial CSI, network)")),
        ],

        "NaiveIntraclusterWMMSE" => [
            ("weighted_logdet_rates_full", Dict(color => "c", :linestyle => "-", :label => "NaiveIntraclusterWMMSE (full CSI, sum)")),
            ("weighted_logdet_rates_partial", Dict(color => "c", :linestyle => "-", :label => "NaiveIntraclusterWMMSE (partial CSI, sum)")),
            ("weighted_logdet_rates_cluster_sdma_full", Dict(color => "c", :linestyle => ":", :label => "NaiveIntraclusterWMMSE (full CSI, cluster)")),
            ("weighted_logdet_rates_network_sdma_full", Dict(color => "c", :linestyle => "--", :label => "NaiveIntraclusterWMMSE (full CSI, network)")),
            ("weighted_logdet_rates_network_sdma_partial", Dict(color => "c", :linestyle => "--", :label => "NaiveIntraclusterWMMSE (partial CSI, network)")),
        ],

        "RobustChen2014_MaxSINR" => [
            ("weighted_logdet_rates_full", Dict(color => "g", :linestyle => "-", :label => "RobustChen2014_MaxSINR (full CSI, sum)")),
            ("weighted_logdet_rates_partial", Dict(color => "g", :linestyle => "-", :label => "RobustChen2014_MaxSINR (partial CSI, sum)")),
            ("weighted_logdet_rates_cluster_sdma_full", Dict(color => "g", :linestyle => ":", :label => "RobustChen2014_MaxSINR (full CSI, cluster)")),
            ("weighted_logdet_rates_network_sdma_full", Dict(color => "g", :linestyle => "--", :label => "RobustChen2014_MaxSINR (full CSI, network)")),
            ("weighted_logdet_rates_network_sdma_partial", Dict(color => "g", :linestyle => "--", :label => "RobustChen2014_MaxSINR (partial CSI, network)")),
        ],

        "NaiveChen2014_MaxSINR" => [
            ("weighted_logdet_rates_full", Dict(color => "y", :linestyle => "-", :label => "NaiveChen2014_MaxSINR (full CSI, sum)")),
            ("weighted_logdet_rates_partial", Dict(color => "y", :linestyle => "-", :label => "NaiveChen2014_MaxSINR (partial CSI, sum)")),
            ("weighted_logdet_rates_cluster_sdma_full", Dict(color => "y", :linestyle => ":", :label => "NaiveChen2014_MaxSINR (full CSI, cluster)")),
            ("weighted_logdet_rates_network_sdma_full", Dict(color => "y", :linestyle => "--", :label => "NaiveChen2014_MaxSINR (full CSI, network)")),
            ("weighted_logdet_rates_network_sdma_partial", Dict(color => "y", :linestyle => "--", :label => "NaiveChen2014_MaxSINR (partial CSI, network)")),
        ],
    )
)

##########################################################################
# Plot it
for file_name in parsed_args["file_names"]
    data = load(file_name)
    processed_results = postprocess(data["raw_precoding_results"], data["simulation_params"], plot_params)
    plot(processed_results, data["simulation_params"], plot_params)
end
