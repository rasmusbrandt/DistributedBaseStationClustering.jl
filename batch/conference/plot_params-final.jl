colours = [
    :ExhaustiveSearchClustering => "#e7298a",
    :CoalitionFormationClustering_Individual => "#1b9e77",
    :Chen2014_kmeans => "#7570b3",
    :RandomClustering => "#66a61e",
    :GrandCoalitionClustering => "#e6ab02",
    :NoClustering => "#d95f02"
]

markers = [
    :ExhaustiveSearchClustering => "*",
    :CoalitionFormationClustering_Individual => "o",
    :Chen2014_kmeans => "p",
    :RandomClustering => "s",
    :GrandCoalitionClustering => "v",
    :NoClustering => "^"
]

labels = [
    :ExhaustiveSearchClustering => "Exhaustive search",
    :CoalitionFormationClustering_Individual => "Coalition formation",
    :Chen2014_kmeans => "k-means clustering [X]",
    :RandomClustering => "Random coalitions",
    :GrandCoalitionClustering => "Grand coalition",
    :NoClustering => "Singleton coalitions"
]

postprocess_params_assignment = [
    "objective" => :sumrate,
    "methods" => [
        "ExhaustiveSearchClustering" => [
            ("utilities",),
            ("no_clusters",),
        ],

        "CoalitionFormationClustering_Individual" => [
            ("utilities",),
            ("no_clusters",),
            ("no_searches",),
        ],

        "Chen2014_kmeans" => [
            ("utilities",),
            ("no_clusters",),
        ],

        "GrandCoalitionClustering" => [
            ("utilities",),
            ("no_clusters",),
        ],

        "RandomClustering" => [
            ("utilities",),
            ("no_clusters",),
        ],

        "NoClustering" => [
            ("utilities",),
            ("no_clusters",),
        ],
    ]
]

postprocess_params_precoding = [
    "objective" => :sumrate,
    "methods" => [
        "ExhaustiveSearchClustering" => [
            ("weighted_logdet_rates_LB",),
        ],

        "CoalitionFormationClustering_Individual" => [
            ("weighted_logdet_rates_LB",),
        ],

        "Chen2014_kmeans" => [
            ("weighted_logdet_rates_LB",),
        ],

        "GrandCoalitionClustering" => [
            ("weighted_logdet_rates_LB",),
        ],

        "RandomClustering" => [
            ("weighted_logdet_rates_LB",),
        ],

        "NoClustering" => [
            ("weighted_logdet_rates_LB",),
        ],
    ]
]
