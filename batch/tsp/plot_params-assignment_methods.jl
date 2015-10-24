using Compat

plot_params_instantaneous_sumrate = Dict(
    "plot_name" => "instantaneous-sumrate",

    "objective" => :sum,

    "figure" => Dict(
        :figsize => (8,5),
        :dpi => 125,
    ),

    "axes" => Dict(
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

        "CoalitionFormationClustering_Attach_IgnoreIAFeasibility" => [
            ("weighted_logdet_rates_full", Dict(:color => "LimeGreen", :linestyle => "--", :label => "CoalitionFormationClustering_Attach (non-IA)")),
        ],

        "CoalitionFormationClustering_AttachOrSupplant_IgnoreIAFeasibility" => [
            ("weighted_logdet_rates_full", Dict(:color => "DarkGreen", :linestyle => "--", :label => "CoalitionFormationClustering_AttachOrSupplant (non-IA)")),
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
    )
)
plot_params_longterm_avg_cluster_size = @compat Dict(
    "plot_name" => "longterm-avg_cluster_size",

    "objective" => :none,

    "figure" => Dict(
        :figsize => (8,5),
        :dpi => 125,
    ),

    "axes" => Dict(
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
plot_params_longterm_num_sum_throughput_calculations = @compat Dict(
    "plot_name" => "longterm-num_sum_throughput_calculations",

    "objective" => :none,

    "figure" => Dict(
        :figsize => (8,5),
        :dpi => 125,
    ),

    "axes" => Dict(
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
