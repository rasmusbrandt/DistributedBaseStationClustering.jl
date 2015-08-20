raw_idx = 1; mean_idx = 2; var_idx = 3

colours_assignment = [
    :BranchAndBoundClustering => "#e7298a",
    :CoalitionFormationClustering_AttachOrSupplant => "#1b9e77",
    :Chen2014_kmeans => "#7570b3",
    :RandomClustering => "#66a61e",
    :GrandCoalitionClustering => "#e6ab02",
    :NoClustering => "#d95f02"
]

linestyles_assignment = [
    :BranchAndBoundClustering => "-",
    :CoalitionFormationClustering_AttachOrSupplant => "-",
    :Chen2014_kmeans => "-",
    :RandomClustering => "-",
    :GrandCoalitionClustering => "-",
    :NoClustering => "-"
]

markers_assignment = [
    :BranchAndBoundClustering => "*",
    :CoalitionFormationClustering_AttachOrSupplant => "o",
    :Chen2014_kmeans => "+",
    :RandomClustering => "s",
    :GrandCoalitionClustering => "v",
    :NoClustering => "^"
]

labels_assignment = [
    :BranchAndBoundClustering => "Global optimum",
    :CoalitionFormationClustering_AttachOrSupplant => "Coalition formation",
    :Chen2014_kmeans => "k-means clustering [x]",
    :RandomClustering => "Random coalitions",
    :GrandCoalitionClustering => "Grand coalition",
    :NoClustering => "Singleton coalitions"
]

postprocess_params_assignment = [
    "objective" => :sumrate,
    "methods" => [
        "BranchAndBoundClustering" => [
            ("throughputs",),
            ("avg_cluster_size",),
        ],

        "CoalitionFormationClustering_AttachOrSupplant" => [
            ("throughputs",),
            ("avg_cluster_size",),
            ("num_searches",),
        ],

        "Chen2014_kmeans" => [
            ("throughputs",),
            ("avg_cluster_size",),
        ],

        "GrandCoalitionClustering" => [
            ("throughputs",),
            ("avg_cluster_size",),
        ],

        "RandomClustering" => [
            ("throughputs",),
            ("avg_cluster_size",),
        ],

        "NoClustering" => [
            ("throughputs",),
            ("avg_cluster_size",),
        ],
    ]
]

postprocess_params_precoding = [
    "objective" => :sumrate,
    "methods" => [
        "BranchAndBoundClustering" => [
            ("weighted_logdet_rates_full",),
        ],

        "CoalitionFormationClustering_AttachOrSupplant" => [
            ("weighted_logdet_rates_full",),
        ],

        "Chen2014_kmeans" => [
            ("weighted_logdet_rates_full",),
        ],

        "GrandCoalitionClustering" => [
            ("weighted_logdet_rates_full",),
        ],

        "RandomClustering" => [
            ("weighted_logdet_rates_full",),
        ],

        "NoClustering" => [
            ("weighted_logdet_rates_full",),
        ],
    ]
]
