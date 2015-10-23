#!/usr/bin/env julia

##########################################################################
# plot_system-SNR-assignment.jl
#
# Plots SNR curves, comparing different cluster assignment methods.
##########################################################################

using DistributedBaseStationClustering, CoordinatedPrecoding

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
    "plot_name" => "instantaneous",

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
            ("weighted_logdet_rates_LB", Dict(:color => "Coral", :linestyle => "", :marker => ".", :label => "ExhaustiveSearchClustering")),
        ],

        "BranchAndBoundClustering" => [
            ("weighted_logdet_rates_LB", Dict(:color => "Coral", :linestyle => "-", :label => "BranchAndBoundClustering")),
        ],


        "CoalitionFormationClustering_Group" => [
            ("weighted_logdet_rates_LB", Dict(:color => "ForestGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Group")),
        ],

        "CoalitionFormationClustering_Individual" => [
            ("weighted_logdet_rates_LB", Dict(:color => "LimeGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Individual")),
        ],


        "GrandCoalitionClustering" => [
            ("weighted_logdet_rates_LB", Dict(:color => "Maroon", :linestyle => "-", :label => "GrandCoalitionClustering")),
        ],

        "GreedyClustering_Single" => [
            ("weighted_logdet_rates_LB", Dict(:color => "DarkOrchid", :linestyle => "-", :label => "GreedyClustering_Single")),
        ],


        "Chen2014_LinearObj_ExhaustiveSearch" => [
            ("weighted_logdet_rates_LB", Dict(:color => "DodgerBlue", :linestyle => "-", :label => "Chen2014_LinearObj_ExhaustiveSearch")),
        ],

        "Chen2014_kmeans" => [
            ("weighted_logdet_rates_LB", Dict(:color => "DodgerBlue", :linestyle => "--", :label => "Chen2014_kmeans")),
        ],

        "Peters2012_Heuristic" => [
            ("weighted_logdet_rates_LB", Dict(:color => "GoldenRod", :linestyle => "-", :label => "Peters2012_Heuristic")),
        ],


        "GreedyClustering_Multiple" => [
            ("weighted_logdet_rates_LB", Dict(:color => "DarkOrchid", :linestyle => "--", :label => "GreedyClustering_Multiple")),
        ],

        "RandomClustering" => [
            ("weighted_logdet_rates_LB", Dict(:color => "Khaki", :linestyle => "-", :label => "RandomClustering")),
        ],

        "NoClustering" => [
            ("weighted_logdet_rates_LB", Dict(:color => "Pink", :linestyle => "-", :label => "NoClustering")),
        ],
    )
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
            ("utilities", Dict(:color => "Coral", :linestyle => "", :marker => ".", :label => "ExhaustiveSearchClustering")),
        ],

        "BranchAndBoundClustering" => [
            ("utilities", Dict(:color => "Coral", :linestyle => "-", :label => "BranchAndBoundClustering")),
        ],


        "CoalitionFormationClustering_Group" => [
            ("utilities", Dict(:color => "ForestGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Group")),
        ],

        "CoalitionFormationClustering_Individual" => [
            ("utilities", Dict(:color => "LimeGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Individual")),
        ],


        "GreedyClustering_Single" => [
            ("utilities", Dict(:color => "DarkOrchid", :linestyle => "-", :label => "GreedyClustering_Single")),
        ],

        "GreedyClustering_Multiple" => [
            ("utilities", Dict(:color => "DarkOrchid", :linestyle => "--", :label => "GreedyClustering_Multiple")),
        ],


        "Chen2014_LinearObj_ExhaustiveSearch" => [
            ("utilities", Dict(:color => "DodgerBlue", :linestyle => "-", :label => "Chen2014_LinearObj_ExhaustiveSearch")),
        ],

        "Chen2014_kmeans" => [
            ("utilities", Dict(:color => "DodgerBlue", :linestyle => "--", :label => "Chen2014_kmeans")),
        ],

        "Peters2012_Heuristic" => [
            ("utilities", Dict(:color => "GoldenRod", :linestyle => "-", :label => "Peters2012_Heuristic")),
        ],


        "GrandCoalitionClustering" => [
            ("utilities", Dict(:color => "Maroon", :linestyle => "-", :label => "GrandCoalitionClustering")),
        ],

        "RandomClustering" => [
            ("utilities", Dict(:color => "Khaki", :linestyle => "-", :label => "RandomClustering")),
        ],

        "NoClustering" => [
            ("utilities", Dict(:color => "Pink", :linestyle => "-", :label => "NoClustering")),
        ],
    )
)
plot_params_longterm_num_utility_calculations = @compat Dict(
    "plot_name" => "longterm-num_utility_calculations",

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
            ("num_utility_calculations", Dict(:color => "Coral", :linestyle => "", :marker => ".", :label => "ExhaustiveSearchClustering")),
        ],

        "BranchAndBoundClustering" => [
            ("num_utility_calculations", Dict(:color => "Coral", :linestyle => "-", :label => "BranchAndBoundClustering")),
        ],


        "CoalitionFormationClustering_Group" => [
            ("num_utility_calculations", Dict(:color => "ForestGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Group")),
        ],

        "CoalitionFormationClustering_Individual" => [
            ("num_utility_calculations", Dict(:color => "LimeGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Individual")),
        ],


        "GreedyClustering_Single" => [
            ("num_utility_calculations", Dict(:color => "DarkOrchid", :linestyle => "-", :label => "GreedyClustering_Single")),
        ],

        "GreedyClustering_Multiple" => [
            ("num_utility_calculations", Dict(:color => "DarkOrchid", :linestyle => "--", :label => "GreedyClustering_Multiple")),
        ],


        "GrandCoalitionClustering" => [
            ("num_utility_calculations", Dict(:color => "Maroon", :linestyle => "-", :label => "GrandCoalitionClustering")),
        ],


        "RandomClustering" => [
            ("num_utility_calculations", Dict(:color => "Khaki", :linestyle => "-", :label => "RandomClustering")),
        ],

        "NoClustering" => [
            ("num_utility_calculations", Dict(:color => "Pink", :linestyle => "-", :label => "NoClustering")),
        ],
    )
)
plot_params_longterm_num_longterm_rate_calculations = @compat Dict(
    "plot_name" => "longterm-num_longterm_rate_calculations",

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
            ("num_longterm_rate_calculations", Dict(:color => "Coral", :linestyle => "", :marker => ".", :label => "ExhaustiveSearchClustering")),
        ],

        "BranchAndBoundClustering" => [
            ("num_longterm_rate_calculations", Dict(:color => "Coral", :linestyle => "-", :label => "BranchAndBoundClustering")),
        ],


        "CoalitionFormationClustering_Group" => [
            ("num_longterm_rate_calculations", Dict(:color => "ForestGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Group")),
        ],

        "CoalitionFormationClustering_Individual" => [
            ("num_longterm_rate_calculations", Dict(:color => "LimeGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Individual")),
        ],


        "GreedyClustering_Single" => [
            ("num_longterm_rate_calculations", Dict(:color => "DarkOrchid", :linestyle => "-", :label => "GreedyClustering_Single")),
        ],

        "GreedyClustering_Multiple" => [
            ("num_longterm_rate_calculations", Dict(:color => "DarkOrchid", :linestyle => "--", :label => "GreedyClustering_Multiple")),
        ],


        "GrandCoalitionClustering" => [
            ("num_longterm_rate_calculations", Dict(:color => "Maroon", :linestyle => "-", :label => "GrandCoalitionClustering")),
        ],

        "RandomClustering" => [
            ("num_longterm_rate_calculations", Dict(:color => "Khaki", :linestyle => "-", :label => "RandomClustering")),
        ],

        "NoClustering" => [
            ("num_longterm_rate_calculations", Dict(:color => "Pink", :linestyle => "-", :label => "NoClustering")),
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
        :ylabel => "Number of utility calculations",
        :yscale => "log",
    ),

    "legend" => Dict(
        :loc => "best",
        :fontsize => 10,
    ),

    "methods" => Dict(
        "CoalitionFormationClustering_Individual" => [
            ("num_searches", Dict(:color => "LimeGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Individual")),
        ],
    )
)
plot_params_longterm_clusters = @compat Dict(
    "plot_name" => "longterm-clusters",

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


        "CoalitionFormationClustering_Group" => [
            ("num_clusters", Dict(:color => "ForestGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Group")),
        ],

        "CoalitionFormationClustering_Individual" => [
            ("num_clusters", Dict(:color => "LimeGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Individual")),
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
    )
)
plot_params_num_searches = @compat Dict(
    "plot_name" => "longterm-clusters",

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
        "CoalitionFormationClustering_Individual" => [
            ("num_searches", Dict(:color => "LimeGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Individual")),
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

    processed_results = postprocess(data["raw_assignment_results"], data["simulation_params"], plot_params_longterm_num_utility_calculations)
    plot(processed_results, data["simulation_params"], plot_params_longterm_num_utility_calculations)

    processed_results = postprocess(data["raw_assignment_results"], data["simulation_params"], plot_params_longterm_num_longterm_rate_calculations)
    plot(processed_results, data["simulation_params"], plot_params_longterm_num_longterm_rate_calculations)

    processed_results = postprocess(data["raw_assignment_results"], data["simulation_params"], plot_params_longterm_num_searches)
    plot(processed_results, data["simulation_params"], plot_params_longterm_num_searches)

    processed_results = postprocess(data["raw_assignment_results"], data["simulation_params"], plot_params_num_searches)
    plot(processed_results, data["simulation_params"], plot_params_num_searches)
end
