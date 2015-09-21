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
    :Chen2014_kmeans => "+",
    :RandomClustering => "s",
    :GrandCoalitionClustering => "v",
    :NoClustering => "^"
]

labels = [
    :ExhaustiveSearchClustering => "Exhaustive search",
    :CoalitionFormationClustering_Individual => "Coalition formation",
    :Chen2014_kmeans => "k-means clustering [4]",
    :RandomClustering => "Random coalitions",
    :GrandCoalitionClustering => "Grand coalition",
    :NoClustering => "Singleton coalitions"
]

postprocess_params_assignment = [
    "objective" => :sumrate,
    "methods" => [
        "ExhaustiveSearchClustering" => [
            ("utilities",),
            ("num_clusters",),
        ],

        "CoalitionFormationClustering_Individual" => [
            ("utilities",),
            ("num_clusters",),
            ("num_searches",),
        ],

        "Chen2014_kmeans" => [
            ("utilities",),
            ("num_clusters",),
        ],

        "GrandCoalitionClustering" => [
            ("utilities",),
            ("num_clusters",),
        ],

        "RandomClustering" => [
            ("utilities",),
            ("num_clusters",),
        ],

        "NoClustering" => [
            ("utilities",),
            ("num_clusters",),
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
