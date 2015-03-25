##########################################################################
# IAClustering
#
# Evaluation environment for the IAClustering project
# https://gitr.sys.kth.se/rabr5411/IAClustering.jl
##########################################################################

module IAClustering

using CoordinatedPrecoding
import Lumberjack

export
    # assignment
    BranchAndBoundClustering,
    Chen2014_ExhaustiveSearch,
    ExhaustiveSearchClustering,
    GrandCoalitionClustering,
    GreedyClustering,
    RandomClustering,
    NoClustering,

    # precoding
    RobustChen2014_MaxSINR, NaiveChen2014_MaxSINR,
    RobustIntraclusterLeakageMinimization, NaiveIntraclusterLeakageMinimization,
    RobustIntraclusterWMMSE, NaiveIntraclusterWMMSE

include("assignment/assignment.jl")

include("precoding/Chen2014_MaxSINR.jl")
include("precoding/IntraclusterLeakageMinimization.jl")
include("precoding/IntraclusterWMMSE.jl")

include("misc/combinations.jl")
include("misc/partitions.jl")
include("misc/feasibility.jl")
include("misc/subsets.jl")
include("misc/utilities.jl")

##########################################################################
# Logging defaults
let
    console = Lumberjack._lumber_mill.timber_trucks["console"]
    Lumberjack.configure(console; mode = "warn")
    file = Lumberjack.add_truck(Lumberjack.LumberjackTruck("default.log", "info"), "default")
end

end
