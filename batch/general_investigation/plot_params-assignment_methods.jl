plot_params_instantaneous_sumrate = [
    "plot_name" => "instantaneous-sumrate",

    "objective" => :sumrate,

    "figure" => [
        :figsize => (8,5),
        :dpi => 125,
    ],

    "axes" => [
        :ylabel => "Sum rate [bits/s/Hz]",
    ],

    "legend" => [
        :loc => "best",
        :fontsize => 10,
    ],

    "methods" => [
        "ExhaustiveSearchClustering" => [
            ("weighted_logdet_rates_full", [ :color => "Coral", :linestyle => "", :marker => ".", :label => "ExhaustiveSearchClustering" ]),
        ],

        "BranchAndBoundClustering" => [
            ("weighted_logdet_rates_full", [ :color => "Coral", :linestyle => "-", :label => "BranchAndBoundClustering" ]),
        ],

        "CoalitionFormationClustering_Attach" => [
            ("weighted_logdet_rates_full", [ :color => "LimeGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Attach" ]),
        ],

        "CoalitionFormationClustering_AttachOrSupplant" => [
            ("weighted_logdet_rates_full", [ :color => "DarkGreen", :linestyle => "-", :label => "CoalitionFormationClustering_AttachOrSupplant" ]),
        ],

        "GrandCoalitionClustering" => [
            ("weighted_logdet_rates_full", [ :color => "Maroon", :linestyle => "-", :label => "GrandCoalitionClustering" ]),
        ],

        "GreedyClustering_Single" => [
            ("weighted_logdet_rates_full", [ :color => "DarkOrchid", :linestyle => "-", :label => "GreedyClustering_Single" ]),
        ],

        "Chen2014_LinearObj_ExhaustiveSearch" => [
            ("weighted_logdet_rates_full", [ :color => "DodgerBlue", :linestyle => "-", :label => "Chen2014_LinearObj_ExhaustiveSearch" ]),
        ],

        "Chen2014_kmeans" => [
            ("weighted_logdet_rates_full", [ :color => "DodgerBlue", :linestyle => "--", :label => "Chen2014_kmeans" ]),
        ],

        "Peters2012_Heuristic" => [
            ("weighted_logdet_rates_full", [ :color => "GoldenRod", :linestyle => "-", :label => "Peters2012_Heuristic" ]),
        ],

        "GreedyClustering_Multiple" => [
            ("weighted_logdet_rates_full", [ :color => "DarkOrchid", :linestyle => "--", :label => "GreedyClustering_Multiple" ]),
        ],

        "RandomClustering" => [
            ("weighted_logdet_rates_full", [ :color => "Khaki", :linestyle => "-", :label => "RandomClustering" ]),
        ],

        "NoClustering" => [
            ("weighted_logdet_rates_full", [ :color => "Pink", :linestyle => "-", :label => "NoClustering" ]),
        ],
    ]
]
plot_params_longterm_sumrate = [
    "plot_name" => "longterm-sumrate",

    "objective" => :sumrate,

    "figure" => [
        :figsize => (8,5),
        :dpi => 125,
    ],

    "axes" => [
        :ylabel => "Sum rate [bits/s/Hz]",
    ],

    "legend" => [
        :loc => "best",
        :fontsize => 10,
    ],

    "methods" => [
        "ExhaustiveSearchClustering" => [
            ("throughputs", [ :color => "Coral", :linestyle => "", :marker => ".", :label => "ExhaustiveSearchClustering" ]),
        ],

        "BranchAndBoundClustering" => [
            ("throughputs", [ :color => "Coral", :linestyle => "-", :label => "BranchAndBoundClustering" ]),
        ],

        "CoalitionFormationClustering_Attach" => [
            ("throughputs", [ :color => "LimeGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Attach" ]),
        ],

        "CoalitionFormationClustering_AttachOrSupplant" => [
            ("throughputs", [ :color => "DarkGreen", :linestyle => "-", :label => "CoalitionFormationClustering_AttachOrSupplant" ]),
        ],

        "GreedyClustering_Single" => [
            ("throughputs", [ :color => "DarkOrchid", :linestyle => "-", :label => "GreedyClustering_Single" ]),
        ],

        "GreedyClustering_Multiple" => [
            ("throughputs", [ :color => "DarkOrchid", :linestyle => "--", :label => "GreedyClustering_Multiple" ]),
        ],

        "Chen2014_LinearObj_ExhaustiveSearch" => [
            ("throughputs", [ :color => "DodgerBlue", :linestyle => "-", :label => "Chen2014_LinearObj_ExhaustiveSearch" ]),
        ],

        "Chen2014_kmeans" => [
            ("throughputs", [ :color => "DodgerBlue", :linestyle => "--", :label => "Chen2014_kmeans" ]),
        ],

        "Peters2012_Heuristic" => [
            ("throughputs", [ :color => "GoldenRod", :linestyle => "-", :label => "Peters2012_Heuristic" ]),
        ],

        "GrandCoalitionClustering" => [
            ("throughputs", [ :color => "Maroon", :linestyle => "-", :label => "GrandCoalitionClustering" ]),
        ],

        "RandomClustering" => [
            ("throughputs", [ :color => "Khaki", :linestyle => "-", :label => "RandomClustering" ]),
        ],

        "NoClustering" => [
            ("throughputs", [ :color => "Pink", :linestyle => "-", :label => "NoClustering" ]),
        ],
    ]
]
plot_params_longterm_avg_cluster_size = [
    "plot_name" => "longterm-avg_cluster_size",

    "objective" => :none,

    "figure" => [
        :figsize => (8,5),
        :dpi => 125,
    ],

    "axes" => [
        :ylabel => "Average cluster size",
    ],

    "legend" => [
        :loc => "best",
        :fontsize => 10,
    ],

    "methods" => [
        "ExhaustiveSearchClustering" => [
            ("avg_cluster_size", [ :color => "Coral", :linestyle => "", :marker => ".", :label => "ExhaustiveSearchClustering" ]),
        ],

        "BranchAndBoundClustering" => [
            ("avg_cluster_size", [ :color => "Coral", :linestyle => "-", :label => "BranchAndBoundClustering" ]),
        ],

        "CoalitionFormationClustering_Attach" => [
            ("avg_cluster_size", [ :color => "LimeGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Attach" ]),
        ],

        "CoalitionFormationClustering_AttachOrSupplant" => [
            ("avg_cluster_size", [ :color => "DarkGreen", :linestyle => "-", :label => "CoalitionFormationClustering_AttachOrSupplant" ]),
        ],

        "GreedyClustering_Single" => [
            ("avg_cluster_size", [ :color => "DarkOrchid", :linestyle => "-", :label => "GreedyClustering_Single" ]),
        ],

        "GreedyClustering_Multiple" => [
            ("avg_cluster_size", [ :color => "DarkOrchid", :linestyle => "--", :label => "GreedyClustering_Multiple" ]),
        ],

        "Chen2014_LinearObj_ExhaustiveSearch" => [
            ("avg_cluster_size", [ :color => "DodgerBlue", :linestyle => "-", :label => "Chen2014_LinearObj_ExhaustiveSearch" ]),
        ],

        "Chen2014_kmeans" => [
            ("avg_cluster_size", [ :color => "DodgerBlue", :linestyle => "--", :label => "Chen2014_kmeans" ]),
        ],

        "GrandCoalitionClustering" => [
            ("avg_cluster_size", [ :color => "Maroon", :linestyle => "-", :label => "GrandCoalitionClustering" ]),
        ],

        "RandomClustering" => [
            ("avg_cluster_size", [ :color => "Khaki", :linestyle => "-", :label => "RandomClustering" ]),
        ],

        "NoClustering" => [
            ("avg_cluster_size", [ :color => "Pink", :linestyle => "-", :label => "NoClustering" ]),
        ],
    ]
]
plot_params_longterm_num_sum_throughput_calculations = [
    "plot_name" => "longterm-num_sum_throughput_calculations",

    "objective" => :none,

    "figure" => [
        :figsize => (8,5),
        :dpi => 125,
    ],

    "axes" => [
        :ylabel => "Number of utility calculations",
        :yscale => "log",
    ],

    "legend" => [
        :loc => "best",
        :fontsize => 10,
    ],

    "methods" => [
        "ExhaustiveSearchClustering" => [
            ("num_sum_throughput_calculations", [ :color => "Coral", :linestyle => "", :marker => ".", :label => "ExhaustiveSearchClustering" ]),
        ],

        "BranchAndBoundClustering" => [
            ("num_sum_throughput_calculations", [ :color => "Coral", :linestyle => "-", :label => "BranchAndBoundClustering" ]),
        ],

        "CoalitionFormationClustering_Attach" => [
            ("num_sum_throughput_calculations", [ :color => "LimeGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Attach" ]),
        ],

        "CoalitionFormationClustering_AttachOrSupplant" => [
            ("num_sum_throughput_calculations", [ :color => "DarkGreen", :linestyle => "-", :label => "CoalitionFormationClustering_AttachOrSupplant" ]),
        ],

        "GreedyClustering_Single" => [
            ("num_sum_throughput_calculations", [ :color => "DarkOrchid", :linestyle => "-", :label => "GreedyClustering_Single" ]),
        ],

        "GreedyClustering_Multiple" => [
            ("num_sum_throughput_calculations", [ :color => "DarkOrchid", :linestyle => "--", :label => "GreedyClustering_Multiple" ]),
        ],

        "GrandCoalitionClustering" => [
            ("num_sum_throughput_calculations", [ :color => "Maroon", :linestyle => "-", :label => "GrandCoalitionClustering" ]),
        ],

        "RandomClustering" => [
            ("num_sum_throughput_calculations", [ :color => "Khaki", :linestyle => "-", :label => "RandomClustering" ]),
        ],

        "NoClustering" => [
            ("num_sum_throughput_calculations", [ :color => "Pink", :linestyle => "-", :label => "NoClustering" ]),
        ],
    ]
]
plot_params_longterm_num_searches = [
    "plot_name" => "longterm-num_searches",

    "objective" => :sumrate,

    "figure" => [
        :figsize => (8,5),
        :dpi => 125,
    ],

    "axes" => [
        :xlabel => "Transmit power [dBm]",
        :ylabel => "Total number of searches",
    ],

    "legend" => [
        :loc => "best",
        :fontsize => 10,
    ],

    "methods" => [
        "CoalitionFormationClustering_Attach" => [
            ("num_searches", [ :color => "LimeGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Attach" ]),
        ],

        "CoalitionFormationClustering_AttachOrSupplant" => [
            ("num_searches", [ :color => "DarkGreen", :linestyle => "-", :label => "CoalitionFormationClustering_AttachOrSupplant" ]),
        ],
    ]
]
