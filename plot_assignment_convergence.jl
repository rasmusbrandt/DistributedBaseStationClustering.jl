#!/usr/bin/env julia

##########################################################################
# plot_assignment_convergence.jl
#
# Plots convergence curves for assignment methods.
##########################################################################

include("src/DistributedBaseStationClustering.jl")
using DistributedBaseStationClustering, CoordinatedPrecoding

##########################################################################
# Load data
#
# Do this before loading other code, otherwise the JLD module might crash!
using Compat, JLD, ArgParse
s = ArgParseSettings()
@add_arg_table s begin
    "file_names"
        help = "file names with results"
        required = true
        nargs = '+'
end
parsed_args = parse_args(s)

##########################################################################
# Plot parameters
plot_params_bounds = @compat Dict(
    "plot_name" => "bounds",

    "objective" => :none,

    "figure" => Dict(
        :figsize => (8,5),
        :dpi => 125,
    ),

    "axes" => Dict(
        :xlabel => "Iterations",
        :ylabel => "Sum rate [bits/s/Hz]",
    ),

    "legend" => Dict(
        :loc => "best",
        :fontsize => 8,
        :ncol => 2,
    Dict(,

    "methods" => Dict(
        "BranchAndBoundClustering" => [
            ("upper_bound_evolution", Dict(:color => "Coral", :linestyle => "--", :label => "BranchAndBoundClustering (UB)")),
            ("lower_bound_evolution", Dict(:color => "Coral", :linestyle => "-", :label => "BranchAndBoundClustering (LB)")),
        ],
    )
)

plot_params_fathom = @compat Dict(
    "plot_name" => "fathom",

    "objective" => :none,

    "figure" => Dict(
        :figsize => (8,5),
        :dpi => 125,
    ),

    "axes" => Dict(
        :xlabel => "Iterations",
        :ylabel => "Fathomed subtree sizes",
    ),

    "legend" => Dict(
        :loc => "best",
        :fontsize => 8,
        :ncol => 2,
    ),

    "methods" => Dict(
        "BranchAndBoundClustering" => [
            ("fathoming_evolution", Dict(:color => "Coral", :linestyle => "--", :label => "BranchAndBoundClustering")),
        ],
    )
)

##########################################################################
# Plot it
for file_name in parsed_args["file_names"]
    data = load(file_name)
    processed_results = postprocess_assignment_convergence(data["raw_results"], data["simulation_params"], plot_params_bounds)
    plot_assignment_convergence(processed_results, data["simulation_params"], plot_params_bounds)
    processed_results = postprocess_assignment_convergence(data["raw_results"], data["simulation_params"], plot_params_fathom)
    plot_assignment_convergence(processed_results, data["simulation_params"], plot_params_fathom)
end
