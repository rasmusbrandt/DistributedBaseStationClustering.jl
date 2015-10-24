using Compat

plot_params = @compat Dict(
    "plot_name" => "",

    "objective" => :sum,

    "figure" => Dict(
        :figsize => (8,5),
        :dpi => 125,
    ),

    "axes" => Dict(
        :ylabel => "Average sum rate [bits/s/Hz]",
    ),

    "legend" => Dict(
        :loc => "best",
        :fontsize => 6,
    ),

    "methods" => Dict(
        "RobustIntraclusterWMMSE" => [
            ("weighted_logdet_rates_full", Dict(:color => "m", :linestyle => "-", :label => "RobustIntraclusterWMMSE (full CSI)")),
            ("weighted_logdet_rates_partial", Dict(:color => "m", :linestyle => "--", :label => "RobustIntraclusterWMMSE (partial CSI)")),
        ],

        "NaiveIntraclusterWMMSE" => [
            ("weighted_logdet_rates_full", Dict(:color => "c", :linestyle => "-", :label => "NaiveIntraclusterWMMSE (full CSI)")),
            ("weighted_logdet_rates_partial", Dict(:color => "c", :linestyle => "--", :label => "NaiveIntraclusterWMMSE (partial CSI)")),
        ],

        "RobustChen2014_MaxSINR" => [
            ("weighted_logdet_rates_full", Dict(:color => "g", :linestyle => "-", :label => "RobustChen2014_MaxSINR (full CSI)")),
            ("weighted_logdet_rates_partial", Dict(:color => "g", :linestyle => "--", :label => "RobustChen2014_MaxSINR (partial CSI)")),
        ],

        "NaiveChen2014_MaxSINR" => [
            ("weighted_logdet_rates_full", Dict(:color => "y", :linestyle => "-", :label => "NaiveChen2014_MaxSINR (full CSI)")),
            ("weighted_logdet_rates_partial", Dict(:color => "y", :linestyle => "--", :label => "NaiveChen2014_MaxSINR (partial CSI)")),
        ],
    )
)
