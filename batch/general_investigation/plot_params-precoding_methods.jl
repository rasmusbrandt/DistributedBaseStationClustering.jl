plot_params = [
    "plot_name" => "",

    "objective" => :sum,

    "figure" => [
        :figsize => (8,5),
        :dpi => 125,
    ],

    "axes" => [
        :ylabel => "Average sum rate [bits/s/Hz]",
    ],

    "legend" => [
        :loc => "best",
        :fontsize => 6,
    ],

    "methods" => [
        "RobustIntraclusterWMMSE" => [
            ("weighted_logdet_rates_full", [ :color => "m", :linestyle => "-", :label => "RobustIntraclusterWMMSE (full CSI)" ]),
            ("weighted_logdet_rates_partial", [ :color => "m", :linestyle => "--", :label => "RobustIntraclusterWMMSE (partial CSI)" ]),
        ],

        "NaiveIntraclusterWMMSE" => [
            ("weighted_logdet_rates_full", [ :color => "c", :linestyle => "-", :label => "NaiveIntraclusterWMMSE (full CSI)" ]),
            ("weighted_logdet_rates_partial", [ :color => "c", :linestyle => "--", :label => "NaiveIntraclusterWMMSE (partial CSI)" ]),
        ],

        "RobustChen2014_MaxSINR" => [
            ("weighted_logdet_rates_full", [ :color => "g", :linestyle => "-", :label => "RobustChen2014_MaxSINR (full CSI)" ]),
            ("weighted_logdet_rates_partial", [ :color => "g", :linestyle => "--", :label => "RobustChen2014_MaxSINR (partial CSI)" ]),
        ],

        "NaiveChen2014_MaxSINR" => [
            ("weighted_logdet_rates_full", [ :color => "y", :linestyle => "-", :label => "NaiveChen2014_MaxSINR (full CSI)" ]),
            ("weighted_logdet_rates_partial", [ :color => "y", :linestyle => "--", :label => "NaiveChen2014_MaxSINR (partial CSI)" ]),
        ],
    ]
]