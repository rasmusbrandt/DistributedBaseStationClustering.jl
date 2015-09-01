plot_params_instantaneous_full_sumrate = [
    "plot_name" => "instantaneous-full",

    "objective" => :sumrate,

    "figure" => [
        :figsize => (8,5),
        :dpi => 125,
    ],

    "axes" => [
        :ylabel => "Average sum rate [bits/s/Hz]",
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


        "CoalitionFormationClustering_Group" => [
            ("weighted_logdet_rates_full", [ :color => "ForestGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Group" ]),
        ],

        "CoalitionFormationClustering_Individual" => [
            ("weighted_logdet_rates_full", [ :color => "LimeGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Individual" ]),
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
plot_params_instantaneous_partial_sumrate = [
    "plot_name" => "instantaneous-partial",

    "objective" => :sumrate,

    "figure" => [
        :figsize => (8,5),
        :dpi => 125,
    ],

    "axes" => [
        :ylabel => "Average sum rate [bits/s/Hz]",
    ],

    "legend" => [
        :loc => "best",
        :fontsize => 10,
    ],

    "methods" => [
        "ExhaustiveSearchClustering" => [
            ("weighted_logdet_rates_partial", [ :color => "Coral", :linestyle => "", :marker => ".", :label => "ExhaustiveSearchClustering" ]),
        ],

        "BranchAndBoundClustering" => [
            ("weighted_logdet_rates_partial", [ :color => "Coral", :linestyle => "-", :label => "BranchAndBoundClustering" ]),
        ],


        "CoalitionFormationClustering_Group" => [
            ("weighted_logdet_rates_partial", [ :color => "ForestGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Group" ]),
        ],

        "CoalitionFormationClustering_Individual" => [
            ("weighted_logdet_rates_partial", [ :color => "LimeGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Individual" ]),
        ],


        "GrandCoalitionClustering" => [
            ("weighted_logdet_rates_partial", [ :color => "Maroon", :linestyle => "-", :label => "GrandCoalitionClustering" ]),
        ],

        "GreedyClustering_Single" => [
            ("weighted_logdet_rates_partial", [ :color => "DarkOrchid", :linestyle => "-", :label => "GreedyClustering_Single" ]),
        ],


        "Chen2014_LinearObj_ExhaustiveSearch" => [
            ("weighted_logdet_rates_partial", [ :color => "DodgerBlue", :linestyle => "-", :label => "Chen2014_LinearObj_ExhaustiveSearch" ]),
        ],

        "Peters2012_Heuristic" => [
            ("weighted_logdet_rates_partial", [ :color => "GoldenRod", :linestyle => "-", :label => "Peters2012_Heuristic" ]),
        ],


        "GreedyClustering_Multiple" => [
            ("weighted_logdet_rates_partial", [ :color => "DarkOrchid", :linestyle => "--", :label => "GreedyClustering_Multiple" ]),
        ],

        "RandomClustering" => [
            ("weighted_logdet_rates_partial", [ :color => "Khaki", :linestyle => "-", :label => "RandomClustering" ]),
        ],

        "NoClustering" => [
            ("weighted_logdet_rates_partial", [ :color => "Pink", :linestyle => "-", :label => "NoClustering" ]),
        ],
    ]
]
plot_params_instantaneous_LB_sumrate = [
    "plot_name" => "instantaneous-LB",

    "objective" => :sumrate,

    "figure" => [
        :figsize => (8,5),
        :dpi => 125,
    ],

    "axes" => [
        :ylabel => "Average sum rate [bits/s/Hz]",
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


        "GrandCoalitionClustering" => [
            ("weighted_logdet_rates_LB", [ :color => "Maroon", :linestyle => "-", :label => "GrandCoalitionClustering" ]),
        ],

        "GreedyClustering_Single" => [
            ("weighted_logdet_rates_LB", [ :color => "DarkOrchid", :linestyle => "-", :label => "GreedyClustering_Single" ]),
        ],


        "Chen2014_LinearObj_ExhaustiveSearch" => [
            ("weighted_logdet_rates_LB", [ :color => "DodgerBlue", :linestyle => "-", :label => "Chen2014_LinearObj_ExhaustiveSearch" ]),
        ],

        "Peters2012_Heuristic" => [
            ("weighted_logdet_rates_LB", [ :color => "GoldenRod", :linestyle => "-", :label => "Peters2012_Heuristic" ]),
        ],


        "GreedyClustering_Multiple" => [
            ("weighted_logdet_rates_LB", [ :color => "DarkOrchid", :linestyle => "--", :label => "GreedyClustering_Multiple" ]),
        ],

        "RandomClustering" => [
            ("weighted_logdet_rates_LB", [ :color => "Khaki", :linestyle => "-", :label => "RandomClustering" ]),
        ],

        "NoClustering" => [
            ("weighted_logdet_rates_LB", [ :color => "Pink", :linestyle => "-", :label => "NoClustering" ]),
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
        :ylabel => "Longterm sum rate (lower bound) [bits/s/Hz]",
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


        "Chen2014_LinearObj_ExhaustiveSearch" => [
            ("utilities", [ :color => "DodgerBlue", :linestyle => "-", :label => "Chen2014_LinearObj_ExhaustiveSearch" ]),
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
    ]
]
plot_params_longterm_num_utility_calculations = [
    "plot_name" => "longterm-num_utility_calculations",

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
            ("num_utility_calculations", [ :color => "Coral", :linestyle => "", :marker => ".", :label => "ExhaustiveSearchClustering" ]),
        ],

        "BranchAndBoundClustering" => [
            ("num_utility_calculations", [ :color => "Coral", :linestyle => "-", :label => "BranchAndBoundClustering" ]),
        ],


        "CoalitionFormationClustering_Group" => [
            ("num_utility_calculations", [ :color => "ForestGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Group" ]),
        ],

        "CoalitionFormationClustering_Individual" => [
            ("num_utility_calculations", [ :color => "LimeGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Individual" ]),
        ],


        "GreedyClustering_Single" => [
            ("num_utility_calculations", [ :color => "DarkOrchid", :linestyle => "-", :label => "GreedyClustering_Single" ]),
        ],

        "GreedyClustering_Multiple" => [
            ("num_utility_calculations", [ :color => "DarkOrchid", :linestyle => "--", :label => "GreedyClustering_Multiple" ]),
        ],


        "GrandCoalitionClustering" => [
            ("num_utility_calculations", [ :color => "Maroon", :linestyle => "-", :label => "GrandCoalitionClustering" ]),
        ],


        "RandomClustering" => [
            ("num_utility_calculations", [ :color => "Khaki", :linestyle => "-", :label => "RandomClustering" ]),
        ],

        "NoClustering" => [
            ("num_utility_calculations", [ :color => "Pink", :linestyle => "-", :label => "NoClustering" ]),
        ],
    ]
]
plot_params_longterm_num_clusters = [
    "plot_name" => "longterm-num_clusters",

    "objective" => :none,

    "figure" => [
        :figsize => (8,5),
        :dpi => 125,
    ],

    "axes" => [
        :ylabel => "Number of clusters",
    ],

    "legend" => [
        :loc => "best",
        :fontsize => 10,
    ],

    "methods" => [
        "ExhaustiveSearchClustering" => [
            ("num_clusters", [ :color => "Coral", :linestyle => "", :marker => ".", :label => "ExhaustiveSearchClustering" ]),
        ],

        "BranchAndBoundClustering" => [
            ("num_clusters", [ :color => "Coral", :linestyle => "-", :label => "BranchAndBoundClustering" ]),
        ],


        "CoalitionFormationClustering_Group" => [
            ("num_clusters", [ :color => "ForestGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Group" ]),
        ],

        "CoalitionFormationClustering_Individual" => [
            ("num_clusters", [ :color => "LimeGreen", :linestyle => "-", :label => "CoalitionFormationClustering_Individual" ]),
        ],


        "GreedyClustering_Single" => [
            ("num_clusters", [ :color => "DarkOrchid", :linestyle => "-", :label => "GreedyClustering_Single" ]),
        ],

        "GreedyClustering_Multiple" => [
            ("num_clusters", [ :color => "DarkOrchid", :linestyle => "--", :label => "GreedyClustering_Multiple" ]),
        ],

        "Chen2014_LinearObj_ExhaustiveSearch" => [
            ("num_clusters", [ :color => "DodgerBlue", :linestyle => "-", :label => "Chen2014_LinearObj_ExhaustiveSearch" ]),
        ],


        "GrandCoalitionClustering" => [
            ("num_clusters", [ :color => "Maroon", :linestyle => "-", :label => "GrandCoalitionClustering" ]),
        ],

        "RandomClustering" => [
            ("num_clusters", [ :color => "Khaki", :linestyle => "-", :label => "RandomClustering" ]),
        ],

        "NoClustering" => [
            ("num_clusters", [ :color => "Pink", :linestyle => "-", :label => "NoClustering" ]),
        ],
    ]
]
