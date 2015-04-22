##########################################################################
# Plot parameters
plot_params_instantaneous_sumrate = [
    "plot_name" => "instantaneous-LB",

    "objective" => :sumrate,

    "figure" => [
        :figsize => (8,5),
        :dpi => 125,
    ],

    "axes" => [
        :xlabel => "Transmit power [dBm]",
        :ylabel => "Sum rate [bits/s/Hz]",
    ],

    "legend" => [
        :loc => "best",
        :fontsize => 10,
    ],

    "methods" => [
        "ExhaustiveSearchClustering" => [
            ("weighted_logdet_rates_LB", [ :color => "Coral", :linestyle => "", :marker => ".", :label => "ExhaustiveSearchClustering" ]),
        ],

        "BranchAndBoundClustering" => [
            ("weighted_logdet_rates_LB", [ :color => "Coral", :linestyle => "-", :label => "BranchAndBoundClustering" ]),
        ],


        "CoalitionFormationClustering_Group" => [
            ("weighted_logdet_rates_LB", [ :color => "ForestGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Group" ]),
        ],

        "CoalitionFormationClustering_Individual" => [
            ("weighted_logdet_rates_LB", [ :color => "LimeGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Individual" ]),
        ],


        "GreedyClustering_Single" => [
            ("weighted_logdet_rates_LB", [ :color => "DarkOrchid", :linestyle => "-", :label => "GreedyClustering_Single" ]),
        ],

        "GreedyClustering_Multiple" => [
            ("weighted_logdet_rates_LB", [ :color => "DarkOrchid", :linestyle => "--", :label => "GreedyClustering_Multiple" ]),
        ],


        "Chen2014_ExhaustiveSearch" => [
            ("weighted_logdet_rates_LB", [ :color => "DodgerBlue", :linestyle => "-", :label => "Chen2014_ExhaustiveSearch" ]),
        ],

        "Peters2012_Heuristic" => [
            ("weighted_logdet_rates_LB", [ :color => "GoldenRod", :linestyle => "-", :label => "Peters2012_Heuristic" ]),
        ],


        "GrandCoalitionClustering" => [
            ("weighted_logdet_rates_LB", [ :color => "Maroon", :linestyle => "-", :label => "GrandCoalitionClustering" ]),
        ],

        "RandomClustering" => [
            ("weighted_logdet_rates_LB", [ :color => "Khaki", :linestyle => "-", :label => "RandomClustering" ]),
        ],

        "NoClustering" => [
            ("weighted_logdet_rates_LB", [ :color => "Pink", :linestyle => "-", :label => "NoClustering" ]),
        ],
    ],
]
plot_params_longterm_sumrate = [
    "plot_name" => "longterm-sumrate",

    "objective" => :sumrate,

    "figure" => [
        :figsize => (8,5),
        :dpi => 125,
    ],

    "axes" => [
        :xlabel => "Transmit power [dBm]",
        :ylabel => "Sum rate [bits/s/Hz]",
    ],

    "legend" => [
        :loc => "best",
        :fontsize => 10,
    ],

    "methods" => [
        "ExhaustiveSearchClustering" => [
            ("utilities", [ :color => "Coral", :linestyle => "", :marker => ".", :label => "ExhaustiveSearchClustering" ]),
        ],

        "BranchAndBoundClustering" => [
            ("utilities", [ :color => "Coral", :linestyle => "-", :label => "BranchAndBoundClustering" ]),
        ],


        "CoalitionFormationClustering_Group" => [
            ("utilities", [ :color => "ForestGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Group" ]),
        ],

        "CoalitionFormationClustering_Individual" => [
            ("utilities", [ :color => "LimeGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Individual" ]),
        ],


        "GreedyClustering_Single" => [
            ("utilities", [ :color => "DarkOrchid", :linestyle => "-", :label => "GreedyClustering_Single" ]),
        ],

        "GreedyClustering_Multiple" => [
            ("utilities", [ :color => "DarkOrchid", :linestyle => "--", :label => "GreedyClustering_Multiple" ]),
        ],


        "Chen2014_ExhaustiveSearch" => [
            ("utilities", [ :color => "DodgerBlue", :linestyle => "-", :label => "Chen2014_ExhaustiveSearch" ]),
        ],

        "Peters2012_Heuristic" => [
            ("utilities", [ :color => "GoldenRod", :linestyle => "-", :label => "Peters2012_Heuristic" ]),
        ],


        "GrandCoalitionClustering" => [
            ("utilities", [ :color => "Maroon", :linestyle => "-", :label => "GrandCoalitionClustering" ]),
        ],

        "RandomClustering" => [
            ("utilities", [ :color => "Khaki", :linestyle => "-", :label => "RandomClustering" ]),
        ],

        "NoClustering" => [
            ("utilities", [ :color => "Pink", :linestyle => "-", :label => "NoClustering" ]),
        ],
    ],
]
plot_params_longterm_no_utility_calculations = [
    "plot_name" => "clustering-no_utility_calculations",

    "objective" => :none,

    "figure" => [
        :figsize => (8,5),
        :dpi => 125,
    ],

    "axes" => [
        :xlabel => "Transmit power [dBm]",
        :ylabel => "Number of utility calculations",
        :yscale => "log",
    ],

    "legend" => [
        :loc => "best",
        :fontsize => 10,
    ],

    "methods" => [
        "ExhaustiveSearchClustering" => [
            ("no_utility_calculations", [ :color => "Coral", :linestyle => "", :marker => ".", :label => "ExhaustiveSearchClustering" ]),
        ],

        "BranchAndBoundClustering" => [
            ("no_utility_calculations", [ :color => "Coral", :linestyle => "-", :label => "BranchAndBoundClustering" ]),
        ],


        "CoalitionFormationClustering_Group" => [
            ("no_utility_calculations", [ :color => "ForestGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Group" ]),
        ],

        "CoalitionFormationClustering_Individual" => [
            ("no_utility_calculations", [ :color => "LimeGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Individual" ]),
        ],


        "GreedyClustering_Single" => [
            ("no_utility_calculations", [ :color => "DarkOrchid", :linestyle => "-", :label => "GreedyClustering_Single" ]),
        ],

        "GreedyClustering_Multiple" => [
            ("no_utility_calculations", [ :color => "DarkOrchid", :linestyle => "--", :label => "GreedyClustering_Multiple" ]),
        ],


        "GrandCoalitionClustering" => [
            ("no_utility_calculations", [ :color => "Maroon", :linestyle => "-", :label => "GrandCoalitionClustering" ]),
        ],

        "RandomClustering" => [
            ("no_utility_calculations", [ :color => "Khaki", :linestyle => "-", :label => "RandomClustering" ]),
        ],

        "NoClustering" => [
            ("no_utility_calculations", [ :color => "Pink", :linestyle => "-", :label => "NoClustering" ]),
        ],
    ],
]
plot_params_longterm_no_longterm_rate_calculations = [
    "plot_name" => "clustering-no_longterm_rate_calculations",

    "objective" => :none,

    "figure" => [
        :figsize => (8,5),
        :dpi => 125,
    ],

    "axes" => [
        :xlabel => "Transmit power [dBm]",
        :ylabel => "Number of utility calculations",
        :yscale => "log",
    ],

    "legend" => [
        :loc => "best",
        :fontsize => 10,
    ],

    "methods" => [
        "ExhaustiveSearchClustering" => [
            ("no_longterm_rate_calculations", [ :color => "Coral", :linestyle => "", :marker => ".", :label => "ExhaustiveSearchClustering" ]),
        ],

        "BranchAndBoundClustering" => [
            ("no_longterm_rate_calculations", [ :color => "Coral", :linestyle => "-", :label => "BranchAndBoundClustering" ]),
        ],


        "CoalitionFormationClustering_Group" => [
            ("no_longterm_rate_calculations", [ :color => "ForestGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Group" ]),
        ],

        "CoalitionFormationClustering_Individual" => [
            ("no_longterm_rate_calculations", [ :color => "LimeGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Individual" ]),
        ],


        "GreedyClustering_Single" => [
            ("no_longterm_rate_calculations", [ :color => "DarkOrchid", :linestyle => "-", :label => "GreedyClustering_Single" ]),
        ],

        "GreedyClustering_Multiple" => [
            ("no_longterm_rate_calculations", [ :color => "DarkOrchid", :linestyle => "--", :label => "GreedyClustering_Multiple" ]),
        ],


        "GrandCoalitionClustering" => [
            ("no_longterm_rate_calculations", [ :color => "Maroon", :linestyle => "-", :label => "GrandCoalitionClustering" ]),
        ],

        "RandomClustering" => [
            ("no_longterm_rate_calculations", [ :color => "Khaki", :linestyle => "-", :label => "RandomClustering" ]),
        ],

        "NoClustering" => [
            ("no_longterm_rate_calculations", [ :color => "Pink", :linestyle => "-", :label => "NoClustering" ]),
        ],
    ],
]
plot_params_longterm_clusters = [
    "plot_name" => "longterm-clusters",

    "objective" => :none,

    "figure" => [
        :figsize => (8,5),
        :dpi => 125,
    ],

    "axes" => [
        :xlabel => "Transmit power [dBm]",
        :ylabel => "Number of clusters",
    ],

    "legend" => [
        :loc => "best",
        :fontsize => 10,
    ],

    "methods" => [
        "ExhaustiveSearchClustering" => [
            ("no_clusters", [ :color => "Coral", :linestyle => "", :marker => ".", :label => "ExhaustiveSearchClustering" ]),
        ],

        "BranchAndBoundClustering" => [
            ("no_clusters", [ :color => "Coral", :linestyle => "-", :label => "BranchAndBoundClustering" ]),
        ],


        "CoalitionFormationClustering_Group" => [
            ("no_clusters", [ :color => "ForestGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Group" ]),
        ],

        "CoalitionFormationClustering_Individual" => [
            ("no_clusters", [ :color => "LimeGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Individual" ]),
        ],


        "GreedyClustering_Single" => [
            ("no_clusters", [ :color => "DarkOrchid", :linestyle => "-", :label => "GreedyClustering_Single" ]),
        ],

        "GreedyClustering_Multiple" => [
            ("no_clusters", [ :color => "DarkOrchid", :linestyle => "--", :label => "GreedyClustering_Multiple" ]),
        ],


        "GrandCoalitionClustering" => [
            ("no_clusters", [ :color => "Maroon", :linestyle => "-", :label => "GrandCoalitionClustering" ]),
        ],

        "RandomClustering" => [
            ("no_clusters", [ :color => "Khaki", :linestyle => "-", :label => "RandomClustering" ]),
        ],

        "NoClustering" => [
            ("no_clusters", [ :color => "Pink", :linestyle => "-", :label => "NoClustering" ]),
        ],
    ],
]
