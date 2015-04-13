plot_params = [
    "plot_name" => "",

    "objective" => :sumrate,

    "figure" => [
        :figsize => (8,5),
        :dpi => 125,
    ],

    "axes" => [
        :xlabel => "Signal-to-noise ratio [dB]",
        :ylabel => "Average sum rate [bits/s/Hz]",
    ],

    "legend" => [
        :loc => "best",
        :fontsize => 6,
    ],

    "methods" => [
        "RobustIntraclusterWMMSE" => [
            ("weighted_logdet_rates", [ :color => "m", :linestyle => "-", :label => "RobustIntraclusterWMMSE" ]),
            ("utilities", [ :color => "m", :linestyle => "--", :label => "RobustIntraclusterWMMSE (utils.)" ]),
        ],

        "NaiveIntraclusterWMMSE" => [
            ("weighted_logdet_rates", [ :color => "m", :linestyle => ":", :label => "NaiveIntraclusterWMMSE" ]),
            # ("utilities", [ :color => "m", :linestyle => "-.", :label => "NaiveIntraclusterWMMSE (utils.)" ]),
        ],


        "RobustIntraclusterLeakageMinimization" => [
            ("weighted_logdet_rates", [ :color => "gray", :linestyle => "-", :label => "RobustIntraclusterLeakageMinimization" ]),
            ("utilities", [ :color => "gray", :linestyle => "--", :label => "RobustIntraclusterLeakageMinimization (utils.)" ]),
        ],

        "NaiveIntraclusterLeakageMinimization" => [
            ("weighted_logdet_rates", [ :color => "gray", :linestyle => ":", :label => "NaiveIntraclusterLeakageMinimization" ]),
            # ("utilities", [ :color => "gray", :linestyle => "-.", :label => "NaiveIntraclusterLeakageMinimization (utils.)" ]),
        ],


        "RobustChen2014_MaxSINR" => [
            ("weighted_logdet_rates", [ :color => "y", :linestyle => "-", :label => "RobustChen2014_MaxSINR" ]),
            ("utilities", [ :color => "y", :linestyle => "--", :label => "RobustChen2014_MaxSINR (utils.)" ]),
        ],

        "NaiveChen2014_MaxSINR" => [
            ("weighted_logdet_rates", [ :color => "y", :linestyle => ":", :label => "NaiveChen2014_MaxSINR" ]),
            # ("utilities", [ :color => "y", :linestyle => "-.", :label => "NaiveChen2014_MaxSINR (utils.)" ]),
        ],


        "Shi2011_WMMSE" => [
            ("weighted_logdet_rates", [ :color => "b", :linestyle => "-", :label => "WMMSE" ]),
        ],

        "Eigenprecoding" => [
            ("intercell_tdma_weighted_logdet_rates", [ :color => "c", :linestyle => "-", :label => "TDMA" ]),
            ("intracell_tdma_weighted_logdet_rates", [ :color => "c", :linestyle => "-.",  :label => "Intracell TDMA" ]),
            ("uncoord_weighted_logdet_rates", [ :color => "k", :linestyle => "-", :label => "Uncoord. transm." ]),
        ],
    ]
]
