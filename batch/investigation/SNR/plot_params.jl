##########################################################################
# Plot parameters
plot_params_instantaneous_sumrate = [
    "plot_name" => "instantaneous-sumrate",

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
            ("weighted_logdet_rates", [ :color => "Coral", :linestyle => "", :marker => ".", :label => "ExhaustiveSearchClustering" ]),
        ],

        "BranchAndBoundClustering" => [
            ("weighted_logdet_rates", [ :color => "Coral", :linestyle => "-", :label => "BranchAndBoundClustering" ]),
        ],


        "CoalitionFormationClustering_Group" => [
            ("weighted_logdet_rates", [ :color => "ForestGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Group" ]),
        ],

        "CoalitionFormationClustering_Individual" => [
            ("weighted_logdet_rates", [ :color => "LimeGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Individual" ]),
        ],


        "Chen2014_ExhaustiveSearch" => [
            ("weighted_logdet_rates", [ :color => "DodgerBlue", :linestyle => "-", :label => "Chen2014_ExhaustiveSearch" ]),
        ],

        "Peters2012_Heuristic" => [
            ("weighted_logdet_rates", [ :color => "GoldenRod", :linestyle => "-", :label => "Peters2012_Heuristic" ]),
        ],


        "GrandCoalitionClustering" => [
            ("weighted_logdet_rates", [ :color => "Maroon", :linestyle => "-", :label => "GrandCoalitionClustering" ]),
        ],

        "GreedyClustering" => [
            ("weighted_logdet_rates", [ :color => "DarkOrchid", :linestyle => "-", :label => "GreedyClustering" ]),
        ],

        "RandomClustering" => [
            ("weighted_logdet_rates", [ :color => "Khaki", :linestyle => "-", :label => "RandomClustering" ]),
        ],

        "NoClustering" => [
            ("weighted_logdet_rates", [ :color => "Pink", :linestyle => "-", :label => "NoClustering" ]),
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


        "Chen2014_ExhaustiveSearch" => [
            ("utilities", [ :color => "DodgerBlue", :linestyle => "-", :label => "Chen2014_ExhaustiveSearch" ]),
        ],

        "Peters2012_Heuristic" => [
            ("utilities", [ :color => "GoldenRod", :linestyle => "-", :label => "Peters2012_Heuristic" ]),
        ],


        "GrandCoalitionClustering" => [
            ("utilities", [ :color => "Maroon", :linestyle => "-", :label => "GrandCoalitionClustering" ]),
        ],

        "GreedyClustering" => [
            ("utilities", [ :color => "DarkOrchid", :linestyle => "-", :label => "GreedyClustering" ]),
        ],

        "RandomClustering" => [
            ("utilities", [ :color => "Khaki", :linestyle => "-", :label => "RandomClustering" ]),
        ],

        "NoClustering" => [
            ("utilities", [ :color => "Pink", :linestyle => "-", :label => "NoClustering" ]),
        ],
    ],
]
plot_params_longterm_iters = [
    "plot_name" => "clustering-iters",

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


        "GrandCoalitionClustering" => [
            ("no_utility_calculations", [ :color => "Maroon", :linestyle => "-", :label => "GrandCoalitionClustering" ]),
        ],

        "GreedyClustering" => [
            ("no_IA_feasibility_checks", [ :color => "DarkOrchid", :linestyle => "-", :label => "GreedyClustering" ]),
        ],

        "RandomClustering" => [
            ("no_utility_calculations", [ :color => "Khaki", :linestyle => "-", :label => "RandomClustering" ]),
        ],

        "NoClustering" => [
            ("no_utility_calculations", [ :color => "Pink", :linestyle => "-", :label => "NoClustering" ]),
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


        "GrandCoalitionClustering" => [
            ("no_clusters", [ :color => "Maroon", :linestyle => "-", :label => "GrandCoalitionClustering" ]),
        ],

        "GreedyClustering" => [
            ("no_IA_feasibility_checks", [ :color => "DarkOrchid", :linestyle => "-", :label => "GreedyClustering" ]),
        ],

        "RandomClustering" => [
            ("no_clusters", [ :color => "Khaki", :linestyle => "-", :label => "RandomClustering" ]),
        ],

        "NoClustering" => [
            ("no_clusters", [ :color => "Pink", :linestyle => "-", :label => "NoClustering" ]),
        ],
    ],
]
