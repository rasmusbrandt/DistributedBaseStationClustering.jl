#!/usr/bin/env julia

##########################################################################
# plot_precoding_convergence-assignment.jl
#
# Plots convergence curves, comparing different cluster assignment methods.
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
        :fontsize => 8,
        :ncol => 2,
    ),

    "methods" => @Compat.Dict(
        "ExhaustiveSearchClustering" => [
            ("weighted_logdet_rates_LB", @Compat.Dict(:color => "Coral", :linestyle => "", :marker => ".", :label => "ExhaustiveSearchClustering")),
        ],

        "BranchAndBoundClustering" => [
            ("weighted_logdet_rates_LB", @Compat.Dict(:color => "Coral", :linestyle => "-", :label => "BranchAndBoundClustering")),
        ],


        "CoalitionFormationClustering_Group" => [
            ("weighted_logdet_rates_LB", @Compat.Dict(:color => "ForestGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Group")),
        ],

        "CoalitionFormationClustering_Individual" => [
            ("weighted_logdet_rates_LB", @Compat.Dict(:color => "LimeGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Individual")),
        ],


        "GreedyClustering_Single" => [
            ("weighted_logdet_rates_LB", @Compat.Dict(:color => "DarkOrchid", :linestyle => "-", :label => "GreedyClustering_Single")),
        ],

        "GreedyClustering_Multiple" => [
            ("weighted_logdet_rates_LB", @Compat.Dict(:color => "DarkOrchid", :linestyle => "--", :label => "GreedyClustering_Multiple")),
        ],


        "Chen2014_LinearObj_ExhaustiveSearch" => [
            ("weighted_logdet_rates_LB", @Compat.Dict(:color => "DodgerBlue", :linestyle => "-", :label => "Chen2014_LinearObj_ExhaustiveSearch")),
        ],


        "GrandCoalitionClustering" => [
            ("weighted_logdet_rates_LB", @Compat.Dict(:color => "Maroon", :linestyle => "-", :label => "GrandCoalitionClustering")),
        ],

        "RandomClustering" => [
            ("weighted_logdet_rates_LB", @Compat.Dict(:color => "Khaki", :linestyle => "-", :label => "RandomClustering")),
        ],

        "NoClustering" => [
            ("weighted_logdet_rates_LB", @Compat.Dict(:color => "Pink", :linestyle => "-", :label => "NoClustering")),
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
