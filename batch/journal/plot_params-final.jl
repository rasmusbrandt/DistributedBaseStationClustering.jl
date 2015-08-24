raw_idx = 1; mean_idx = 2; var_idx = 3

# 8-class Set1
colours = [
    :red => "#e41a1c",
    :blue => "#377eb8",
    :green => "#4daf4a",
    :purple => "#984ea3",
    :orange => "#ff7f00",
    :yellow => "#ffff33",
    :brown => "#a65628",
    :pink => "#f781bf",
]

colours_assignment = [
    :BranchAndBoundClustering => colours[:red],
    :CoalitionFormationClustering_AttachOrSupplant => colours[:blue],
    :CoalitionFormationClustering_Attach => colours[:green],
    :Peters2012_Heuristic => colours[:purple],
    :Chen2014_kmeans => colours[:pink],
    :RandomClustering => colours[:yellow],
    :GrandCoalitionClustering => colours[:brown],
    :NoClustering => colours[:orange],
]

linestyles_assignment = [
    :BranchAndBoundClustering => "-",
    :CoalitionFormationClustering_AttachOrSupplant => "-",
    :CoalitionFormationClustering_Attach => "-",
    :Chen2014_kmeans => "-",
    :Peters2012_Heuristic => "-",
    :RandomClustering => "-",
    :GrandCoalitionClustering => "-",
    :NoClustering => "-"
]

markers_assignment = [
    :BranchAndBoundClustering => "*",
    :CoalitionFormationClustering_AttachOrSupplant => "o",
    :CoalitionFormationClustering_Attach => "o",
    :Chen2014_kmeans => "+",
    :Peters2012_Heuristic => "+",
    :RandomClustering => "s",
    :GrandCoalitionClustering => "v",
    :NoClustering => "^"
]

labels_assignment = [
    :BranchAndBoundClustering => "Global optimum",
    :CoalitionFormationClustering_AttachOrSupplant => "Coalition formation",
    :CoalitionFormationClustering_Attach => "Coalition formation (attach only)",
    :Chen2014_kmeans => "k-means clustering [x]",
    :Peters2012_Heuristic => "Grouping heuristic [x]",
    :RandomClustering => "Random coalitions",
    :GrandCoalitionClustering => "Grand coalition",
    :NoClustering => "Singleton coalitions"
]

postprocess_params_assignment = [
    "objective" => :sumrate,
    "methods" => [
        "BranchAndBoundClustering" => [
            ("throughputs",),
            ("throughputs_cluster_sdma",),
            ("throughputs_network_sdma",),
            ("avg_cluster_size",),
        ],

        "CoalitionFormationClustering_AttachOrSupplant" => [
            ("throughputs",),
            ("throughputs_cluster_sdma",),
            ("throughputs_network_sdma",),
            ("avg_cluster_size",),
            ("num_searches",),
        ],

        "CoalitionFormationClustering_Attach" => [
            ("throughputs",),
            ("throughputs_cluster_sdma",),
            ("throughputs_network_sdma",),
            ("avg_cluster_size",),
            ("num_searches",),
        ],

        "Chen2014_kmeans" => [
            ("throughputs",),
            ("throughputs_cluster_sdma",),
            ("throughputs_network_sdma",),
            ("avg_cluster_size",),
        ],

        "Peters2012_Heuristic" => [
            ("throughputs",),
            ("throughputs_cluster_sdma",),
            ("throughputs_network_sdma",),
            ("avg_cluster_size",),
        ],

        "GrandCoalitionClustering" => [
            ("throughputs",),
            ("throughputs_cluster_sdma",),
            ("throughputs_network_sdma",),
            ("avg_cluster_size",),
        ],

        "RandomClustering" => [
            ("throughputs",),
            ("throughputs_cluster_sdma",),
            ("throughputs_network_sdma",),
            ("avg_cluster_size",),
        ],

        "NoClustering" => [
            ("throughputs",),
            ("throughputs_cluster_sdma",),
            ("throughputs_network_sdma",),
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

        "CoalitionFormationClustering_Attach" => [
            ("weighted_logdet_rates_full",),
        ],

        "Chen2014_kmeans" => [
            ("weighted_logdet_rates_full",),
        ],

        "Peters2012_Heuristic" => [
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
