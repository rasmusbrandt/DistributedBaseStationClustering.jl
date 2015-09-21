#!/usr/bin/env julia

##########################################################################
# plot_precoding_convergence-precoding.jl
#
# Plots convergence curves, comparing different precoding methods.
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

    "objective" => :sumrate,

    "figure" => @Compat.Dict(
        :figsize => (8,5),
        :dpi => 125,
    ),

    "axes" => @Compat.Dict(
        :xlabel => "Iterations",
        :ylabel => "Sum rate [bits/s/Hz]",
        :ylim => [0, 70],
    ),

    "legend" => @Compat.Dict(
        :loc => "best",
        :fontsize => 6,
        :ncol => 4,
    ),

    "methods" => @Compat.Dict(
        "RobustIntraclusterWMMSE" => [
            ("weighted_logdet_rates_full", @Compat.Dict(:color => "m", :linestyle => "-", :label => "RobustIntraclusterWMMSE (full CSI)")),
            ("weighted_logdet_rates_partial", @Compat.Dict(:color => "m", :linestyle => "--", :label => "RobustIntraclusterWMMSE (partial CSI)")),
            ("weighted_logdet_rates_LB", @Compat.Dict(:color => "m", :linestyle => ":", :label => "RobustIntraclusterWMMSE (lower bound)")),
        ],

        "NaiveIntraclusterWMMSE" => [
            ("weighted_logdet_rates_full", @Compat.Dict(:color => "c", :linestyle => "-", :label => "NaiveIntraclusterWMMSE (full CSI)")),
            ("weighted_logdet_rates_partial", @Compat.Dict(:color => "c", :linestyle => "--", :label => "NaiveIntraclusterWMMSE (partial CSI)")),
            ("weighted_logdet_rates_LB", @Compat.Dict(:color => "c", :linestyle => ":", :label => "NaiveIntraclusterWMMSE (lower bound)")),
        ],


        "RobustIntraclusterLeakageMinimization" => [
            ("weighted_logdet_rates_full", @Compat.Dict(:color => "DarkGray", :linestyle => "-", :label => "RobustIntraclusterLeakageMinimization (full CSI)")),
            ("weighted_logdet_rates_partial", @Compat.Dict(:color => "DarkGray", :linestyle => "--", :label => "RobustIntraclusterLeakageMinimization (partial CSI)")),
            ("weighted_logdet_rates_LB", @Compat.Dict(:color => "DarkGray", :linestyle => ":", :label => "RobustIntraclusterLeakageMinimization (lower bound)")),
        ],

        "NaiveIntraclusterLeakageMinimization" => [
            ("weighted_logdet_rates_full", @Compat.Dict(:color => "gray", :linestyle => "-", :label => "NaiveIntraclusterLeakageMinimization (full CSI)")),
            ("weighted_logdet_rates_partial", @Compat.Dict(:color => "gray", :linestyle => "--", :label => "NaiveIntraclusterLeakageMinimization (partial CSI)")),
            ("weighted_logdet_rates_LB", @Compat.Dict(:color => "gray", :linestyle => ":", :label => "NaiveIntraclusterLeakageMinimization (lower bound)")),
        ],


        "RobustChen2014_MaxSINR" => [
            ("weighted_logdet_rates_full", @Compat.Dict(:color => "g", :linestyle => "-", :label => "RobustChen2014_MaxSINR (full CSI)")),
            ("weighted_logdet_rates_partial", @Compat.Dict(:color => "g", :linestyle => "--", :label => "RobustChen2014_MaxSINR (partial CSI)")),
            ("weighted_logdet_rates_LB", @Compat.Dict(:color => "g", :linestyle => ":", :label => "RobustChen2014_MaxSINR (lower bound)")),
        ],

        "NaiveChen2014_MaxSINR" => [
            ("weighted_logdet_rates_full", @Compat.Dict(:color => "y", :linestyle => "-", :label => "NaiveChen2014_MaxSINR (full CSI)")),
            ("weighted_logdet_rates_partial", @Compat.Dict(:color => "y", :linestyle => "--", :label => "NaiveChen2014_MaxSINR (partial CSI)")),
            ("weighted_logdet_rates_LB", @Compat.Dict(:color => "y", :linestyle => ":", :label => "NaiveChen2014_MaxSINR (lower bound)")),
        ],


        "Shi2011_WMMSE" => [
            ("weighted_logdet_rates", @Compat.Dict(:color => "b", :linestyle => "-", :label => "WMMSE (full CSI)")),
        ],

        "Eigenprecoding" => [
            ("intercell_tdma_weighted_logdet_rates", @Compat.Dict(:color => "k", :linestyle => "-", :label => "TDMA")),
            ("intracell_tdma_weighted_logdet_rates", @Compat.Dict(:color => "k", :linestyle => "--",  :label => "Intracell TDMA")),
            ("uncoord_weighted_logdet_rates", @Compat.Dict(:color => "k", :linestyle => ":", :label => "Uncoord. transm.")),
        ],
    )
)

##########################################################################
# Plot it
for file_name in parsed_args["file_names"]
    data = load(file_name)
    processed_results = postprocess_precoding_convergence(data["raw_results"], data["simulation_params"], plot_params)
    plot_precoding_convergence(processed_results, data["simulation_params"], plot_params)
end
