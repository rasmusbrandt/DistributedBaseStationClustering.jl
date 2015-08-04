##########################################################################
# IAClustering
#
# Evaluation environment for the IAClustering project
# https://gitr.sys.kth.se/rabr5411/IAClustering.jl
##########################################################################

module IAClustering

using CoordinatedPrecoding
import Lumberjack, Clustering, Polynomials

# Since my covariance matrices are small (2x2, 4x4, etc.), there is no
# gain from using multithreaded BLAS!
blas_set_num_threads(1)

export
    # assignment
    BranchAndBoundClustering,
    Chen2014_LinearObj_ExhaustiveSearch,
    Chen2014_kmeans,
    CoalitionFormationClustering_Group,
    CoalitionFormationClustering_Individual,
    CoalitionFormationClustering_Swap,
    ExhaustiveSearchClustering,
    GrandCoalitionClustering,
    GreedyClustering_Single,
    GreedyClustering_Multiple,
    RandomClustering,
    NoClustering,
    Peters2012_Heuristic,

    # precoding
    RobustChen2014_MaxSINR, NaiveChen2014_MaxSINR,
    RobustIntraclusterLeakageMinimization, NaiveIntraclusterLeakageMinimization,
    RobustIntraclusterWMMSE, NaiveIntraclusterWMMSE,
    NoPrecoding

include("misc/combinations.jl")
include("misc/expint.jl")
include("misc/partitions.jl")
include("misc/feasibility.jl")
include("misc/subsets.jl")
include("misc/utilities.jl")

include("assignment/assignment.jl")
include("assignment/BranchAndBoundClustering.jl")
include("assignment/Chen2014.jl")
include("assignment/CoalitionFormationClustering.jl")
include("assignment/ExhaustiveSearchClustering.jl")
include("assignment/GrandCoalitionClustering.jl")
include("assignment/GreedyClustering.jl")
include("assignment/NoClustering.jl")
include("assignment/Peters2012_Heuristic.jl")
include("assignment/RandomClustering.jl")

include("precoding/precoding.jl")
include("precoding/Chen2014_MaxSINR.jl")
include("precoding/IntraclusterLeakageMinimization.jl")
include("precoding/IntraclusterWMMSE.jl")
include("precoding/NoPrecoding.jl")

##########################################################################
# Logging defaults
let
    console = Lumberjack._lumber_mill.timber_trucks["console"]
    Lumberjack.configure(console; mode = "warn")
    file = Lumberjack.add_truck(Lumberjack.LumberjackTruck("default.log", "info"), "default")
end

end
