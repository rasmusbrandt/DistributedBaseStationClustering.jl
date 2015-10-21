#!/usr/bin/env julia

##########################################################################
# plot_system-SNR-assignment.jl
#
# Plots SNR curves, comparing different cluster assignment methods.
##########################################################################

using CoordinatedPrecoding, DistributedBaseStationClustering

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
plot_params_instantaneous_sumrate = @compat Dict(
    "plot_name" => "instantaneous-sumrate",

    "objective" => :sum,

    "figure" => Dict(
        :figsize => (8,5),
        :dpi => 125,
    ),

    "axes" => Dict(
        :xlabel => "SNR [dB]",
        :ylabel => "Sum rate [bits/s/Hz]",
    ),

    "legend" => Dict(
        :loc => "best",
        :fontsize => 10,
    ),

    "methods" => Dict(
        "ExhaustiveSearchClustering" => [
            ("weighted_logdet_rates_full", Dict(:color => "Coral", :linestyle => "", :marker => ".", :label => "ExhaustiveSearchClustering")),
        ],

        "BranchAndBoundClustering" => [
            ("weighted_logdet_rates_full", Dict(:color => "Coral", :linestyle => "-", :label => "BranchAndBoundClustering")),
        ],

        "CoalitionFormationClustering_Attach" => [
            ("weighted_logdet_rates_full", Dict(:color => "LimeGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Attach")),
        ],

        "CoalitionFormationClustering_AttachOrSupplant" => [
            ("weighted_logdet_rates_full", Dict(:color => "DarkGreen", :linestyle => "-", :label => "CoalitionFormationClustering_AttachOrSupplant")),
        ],

        "GrandCoalitionClustering" => [
            ("weighted_logdet_rates_full", Dict(:color => "Maroon", :linestyle => "-", :label => "GrandCoalitionClustering")),
        ],

        "GreedyClustering_Single" => [
            ("weighted_logdet_rates_full", Dict(:color => "DarkOrchid", :linestyle => "-", :label => "GreedyClustering_Single")),
        ],

        "Chen2014_LinearObj_ExhaustiveSearch" => [
            ("weighted_logdet_rates_full", Dict(:color => "DodgerBlue", :linestyle => "-", :label => "Chen2014_LinearObj_ExhaustiveSearch")),
        ],

        "Chen2014_kmeans" => [
            ("weighted_logdet_rates_full", Dict(:color => "DodgerBlue", :linestyle => "--", :label => "Chen2014_kmeans")),
        ],

        "Peters2012_Heuristic" => [
            ("weighted_logdet_rates_full", Dict(:color => "GoldenRod", :linestyle => "-", :label => "Peters2012_Heuristic")),
        ],

        "GreedyClustering_Multiple" => [
            ("weighted_logdet_rates_full", Dict(:color => "DarkOrchid", :linestyle => "--", :label => "GreedyClustering_Multiple")),
        ],

        "RandomClustering" => [
            ("weighted_logdet_rates_full", Dict(:color => "Khaki", :linestyle => "-", :label => "RandomClustering")),
        ],

        "NoClustering" => [
            ("weighted_logdet_rates_full", Dict(:color => "Pink", :linestyle => "-", :label => "NoClustering")),
        ],
    ),
)
plot_params_longterm_sumrate = @compat Dict(
    "plot_name" => "longterm-sumrate",

    "objective" => :sum,

    "figure" => Dict(
        :figsize => (8,5),
        :dpi => 125,
    ),

    "axes" => Dict(
        :xlabel => "Transmit power [dBm]",
        :ylabel => "Sum rate [bits/s/Hz]",
    ),

    "legend" => Dict(
        :loc => "best",
        :fontsize => 10,
    ),

    "methods" => Dict(
        "ExhaustiveSearchClustering" => [
            ("throughputs", Dict(:color => "Coral", :linestyle => "", :marker => ".", :label => "ExhaustiveSearchClustering")),
        ],

        "BranchAndBoundClustering" => [
            ("throughputs", Dict(:color => "Coral", :linestyle => "-", :label => "BranchAndBoundClustering")),
        ],

        "CoalitionFormationClustering_Attach" => [
            ("throughputs", Dict(:color => "LimeGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Attach")),
        ],

        "CoalitionFormationClustering_AttachOrSupplant" => [
            ("throughputs", Dict(:color => "DarkGreen", :linestyle => "-", :label => "CoalitionFormationClustering_AttachOrSupplant")),
        ],

        "GreedyClustering_Single" => [
            ("throughputs", Dict(:color => "DarkOrchid", :linestyle => "-", :label => "GreedyClustering_Single")),
        ],

        "GreedyClustering_Multiple" => [
            ("throughputs", Dict(:color => "DarkOrchid", :linestyle => "--", :label => "GreedyClustering_Multiple")),
        ],

        "Chen2014_LinearObj_ExhaustiveSearch" => [
            ("throughputs", Dict(:color => "DodgerBlue", :linestyle => "-", :label => "Chen2014_LinearObj_ExhaustiveSearch")),
        ],

        "Chen2014_kmeans" => [
            ("throughputs", Dict(:color => "DodgerBlue", :linestyle => "--", :label => "Chen2014_kmeans")),
        ],

        "Peters2012_Heuristic" => [
            ("throughputs", Dict(:color => "GoldenRod", :linestyle => "-", :label => "Peters2012_Heuristic")),
        ],

        "GrandCoalitionClustering" => [
            ("throughputs", Dict(:color => "Maroon", :linestyle => "-", :label => "GrandCoalitionClustering")),
        ],

        "RandomClustering" => [
            ("throughputs", Dict(:color => "Khaki", :linestyle => "-", :label => "RandomClustering")),
        ],

        "NoClustering" => [
            ("throughputs", Dict(:color => "Pink", :linestyle => "-", :label => "NoClustering")),
        ],
    ),
)
plot_params_longterm_num_sum_throughput_calculations = @compat Dict(
    "plot_name" => "longterm-num_sum_throughput_calculations",

    "objective" => :none,

    "figure" => Dict(
        :figsize => (8,5),
        :dpi => 125,
    ),

    "axes" => Dict(
        :xlabel => "Transmit power [dBm]",
        :ylabel => "Number of utility calculations",
        :yscale => "log",
    ),

    "legend" => Dict(
        :loc => "best",
        :fontsize => 10,
    ),

    "methods" => Dict(
        "ExhaustiveSearchClustering" => [
            ("num_sum_throughput_calculations", Dict(:color => "Coral", :linestyle => "", :marker => ".", :label => "ExhaustiveSearchClustering")),
        ],

        "BranchAndBoundClustering" => [
            ("num_sum_throughput_calculations", Dict(:color => "Coral", :linestyle => "-", :label => "BranchAndBoundClustering")),
        ],

        "CoalitionFormationClustering_Attach" => [
            ("num_sum_throughput_calculations", Dict(:color => "LimeGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Attach")),
        ],

        "CoalitionFormationClustering_AttachOrSupplant" => [
            ("num_sum_throughput_calculations", Dict(:color => "DarkGreen", :linestyle => "-", :label => "CoalitionFormationClustering_AttachOrSupplant")),
        ],

        "GreedyClustering_Single" => [
            ("num_sum_throughput_calculations", Dict(:color => "DarkOrchid", :linestyle => "-", :label => "GreedyClustering_Single")),
        ],

        "GreedyClustering_Multiple" => [
            ("num_sum_throughput_calculations", Dict(:color => "DarkOrchid", :linestyle => "--", :label => "GreedyClustering_Multiple")),
        ],

        "GrandCoalitionClustering" => [
            ("num_sum_throughput_calculations", Dict(:color => "Maroon", :linestyle => "-", :label => "GrandCoalitionClustering")),
        ],

        "RandomClustering" => [
            ("num_sum_throughput_calculations", Dict(:color => "Khaki", :linestyle => "-", :label => "RandomClustering")),
        ],

        "NoClustering" => [
            ("num_sum_throughput_calculations", Dict(:color => "Pink", :linestyle => "-", :label => "NoClustering")),
        ],
    ),
)
plot_params_longterm_num_clusters = @compat Dict(
    "plot_name" => "longterm-num_clusters",

    "objective" => :none,

    "figure" => Dict(
        :figsize => (8,5),
        :dpi => 125,
    ),

    "axes" => Dict(
        :xlabel => "Transmit power [dBm]",
        :ylabel => "Number of clusters",
    ),

    "legend" => Dict(
        :loc => "best",
        :fontsize => 10,
    ),

    "methods" => Dict(
        "ExhaustiveSearchClustering" => [
            ("num_clusters", Dict(:color => "Coral", :linestyle => "", :marker => ".", :label => "ExhaustiveSearchClustering")),
        ],

        "BranchAndBoundClustering" => [
            ("num_clusters", Dict(:color => "Coral", :linestyle => "-", :label => "BranchAndBoundClustering")),
        ],

        "CoalitionFormationClustering_Attach" => [
            ("num_clusters", Dict(:color => "LimeGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Attach")),
        ],

        "CoalitionFormationClustering_AttachOrSupplant" => [
            ("num_clusters", Dict(:color => "DarkGreen", :linestyle => "-", :label => "CoalitionFormationClustering_AttachOrSupplant")),
        ],

        "GreedyClustering_Single" => [
            ("num_clusters", Dict(:color => "DarkOrchid", :linestyle => "-", :label => "GreedyClustering_Single")),
        ],

        "GreedyClustering_Multiple" => [
            ("num_clusters", Dict(:color => "DarkOrchid", :linestyle => "--", :label => "GreedyClustering_Multiple")),
        ],

        "Chen2014_LinearObj_ExhaustiveSearch" => [
            ("num_clusters", Dict(:color => "DodgerBlue", :linestyle => "-", :label => "Chen2014_LinearObj_ExhaustiveSearch")),
        ],

        "Chen2014_kmeans" => [
            ("num_clusters", Dict(:color => "DodgerBlue", :linestyle => "--", :label => "Chen2014_kmeans")),
        ],

        "GrandCoalitionClustering" => [
            ("num_clusters", Dict(:color => "Maroon", :linestyle => "-", :label => "GrandCoalitionClustering")),
        ],

        "RandomClustering" => [
            ("num_clusters", Dict(:color => "Khaki", :linestyle => "-", :label => "RandomClustering")),
        ],

        "NoClustering" => [
            ("num_clusters", Dict(:color => "Pink", :linestyle => "-", :label => "NoClustering")),
        ],
    ),
)
plot_params_longterm_avg_cluster_size = @compat Dict(
    "plot_name" => "longterm-avg_cluster_size",

    "objective" => :none,

    "figure" => Dict(
        :figsize => (8,5),
        :dpi => 125,
    ),

    "axes" => Dict(
        :xlabel => "Transmit power [dBm]",
        :ylabel => "Average cluster size",
    ),

    "legend" => Dict(
        :loc => "best",
        :fontsize => 10,
    ),

    "methods" => Dict(
        "ExhaustiveSearchClustering" => [
            ("avg_cluster_size", Dict(:color => "Coral", :linestyle => "", :marker => ".", :label => "ExhaustiveSearchClustering")),
        ],

        "BranchAndBoundClustering" => [
            ("avg_cluster_size", Dict(:color => "Coral", :linestyle => "-", :label => "BranchAndBoundClustering")),
        ],

        "CoalitionFormationClustering_Attach" => [
            ("avg_cluster_size", Dict(:color => "LimeGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Attach")),
        ],

        "CoalitionFormationClustering_AttachOrSupplant" => [
            ("avg_cluster_size", Dict(:color => "DarkGreen", :linestyle => "-", :label => "CoalitionFormationClustering_AttachOrSupplant")),
        ],

        "GreedyClustering_Single" => [
            ("avg_cluster_size", Dict(:color => "DarkOrchid", :linestyle => "-", :label => "GreedyClustering_Single")),
        ],

        "GreedyClustering_Multiple" => [
            ("avg_cluster_size", Dict(:color => "DarkOrchid", :linestyle => "--", :label => "GreedyClustering_Multiple")),
        ],

        "Chen2014_LinearObj_ExhaustiveSearch" => [
            ("avg_cluster_size", Dict(:color => "DodgerBlue", :linestyle => "-", :label => "Chen2014_LinearObj_ExhaustiveSearch")),
        ],

        "Chen2014_kmeans" => [
            ("avg_cluster_size", Dict(:color => "DodgerBlue", :linestyle => "--", :label => "Chen2014_kmeans")),
        ],

        "GrandCoalitionClustering" => [
            ("avg_cluster_size", Dict(:color => "Maroon", :linestyle => "-", :label => "GrandCoalitionClustering")),
        ],

        "RandomClustering" => [
            ("avg_cluster_size", Dict(:color => "Khaki", :linestyle => "-", :label => "RandomClustering")),
        ],

        "NoClustering" => [
            ("avg_cluster_size", Dict(:color => "Pink", :linestyle => "-", :label => "NoClustering")),
        ],
    )
)
plot_params_longterm_num_searches = @compat Dict(
    "plot_name" => "longterm-num_searches",

    "objective" => :sum,

    "figure" => Dict(
        :figsize => (8,5),
        :dpi => 125,
    ),

    "axes" => Dict(
        :xlabel => "Transmit power [dBm]",
        :ylabel => "Total number of searches",
    ),

    "legend" => Dict(
        :loc => "best",
        :fontsize => 10,
    ),

    "methods" => Dict(
        "CoalitionFormationClustering_Attach" => [
            ("num_searches", Dict(:color => "LimeGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Attach")),
        ],

        "CoalitionFormationClustering_AttachOrSupplant" => [
            ("num_searches", Dict(:color => "DarkGreen", :linestyle => "-", :label => "CoalitionFormationClustering_AttachOrSupplant")),
        ],
    )
)

##########################################################################
# Plot it
for file_name in parsed_args["file_names"]
    data = load(file_name)

    processed_results = postprocess(data["raw_precoding_results"], data["simulation_params"], plot_params_instantaneous_sumrate)
    plot(processed_results, data["simulation_params"], plot_params_instantaneous_sumrate)

    processed_results = postprocess(data["raw_assignment_results"], data["simulation_params"], plot_params_longterm_sumrate)
    plot(processed_results, data["simulation_params"], plot_params_longterm_sumrate)

    processed_results = postprocess(data["raw_assignment_results"], data["simulation_params"], plot_params_longterm_num_sum_throughput_calculations)
    plot(processed_results, data["simulation_params"], plot_params_longterm_num_sum_throughput_calculations)

    processed_results = postprocess(data["raw_assignment_results"], data["simulation_params"], plot_params_longterm_num_clusters)
    plot(processed_results, data["simulation_params"], plot_params_longterm_num_clusters)

    processed_results = postprocess(data["raw_assignment_results"], data["simulation_params"], plot_params_longterm_avg_cluster_size)
    plot(processed_results, data["simulation_params"], plot_params_longterm_avg_cluster_size)

    processed_results = postprocess(data["raw_assignment_results"], data["simulation_params"], plot_params_longterm_num_searches)
    plot(processed_results, data["simulation_params"], plot_params_longterm_num_searches)
end
