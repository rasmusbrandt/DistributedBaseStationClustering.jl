##########################################################################
# DistributedBSClustering
#
# Evaluation environment for the DistributedBSClustering project
# https://gitr.sys.kth.se/rabr5411/DistributedBSClustering.jl
##########################################################################

module DistributedBSClustering

using CoordinatedPrecoding
import Lumberjack, PyCall, Clustering
using Compat

PyCall.@pyimport scipy.special as scipy_special

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
    RobustIntraclusterWMMSE, NaiveIntraclusterWMMSE

include("misc/combinations.jl")
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

##########################################################################
# Logging defaults
let
    console = Lumberjack._lumber_mill.timber_trucks["console"]
    Lumberjack.configure(console; mode = "warn")
    file = Lumberjack.add_truck(Lumberjack.LumberjackTruck("default.log", "info"), "default")
end

##########################################################################
# Hermitian methods used in precoding methods

if VERSION >= v"0.4-pre"
    macro hermdata(A)
        return :($A.data)
    end
else
    macro hermdata(A)
        return :($A.S)
    end
end

import Base: +, -, .*
+(A::Hermitian{Complex128}, B::Hermitian{Complex128}) = Hermitian(hermdata(A) + hermdata(B))
+(B::Matrix{Float64}, A::Hermitian{Complex128}) = +(A, B)
+(A::Hermitian{Complex128}, B::Matrix{Float64}) = hermdata(A) + B
+(A::Hermitian{Complex128}, B::Matrix{Complex128}) = hermdata(A) + B

-(A::Hermitian{Complex128}, B::Matrix{Complex128}) = hermdata(A) - B
-(B::Array{Complex128, 2}, A::Hermitian{Complex128}) = -(A, B)
-(A::Hermitian{Complex128}, B::Matrix{Float64}) = hermdata(A) - B

.*(a::Float64, B::Hermitian{Complex128}) = Hermitian(a.*hermdata(B))

import Base.logdet, Base.diag
logdet(A::Hermitian{Complex128}) = logdet(hermdata(A))
diag(A::Hermitian{Complex128}) = diag(hermdata(A))

end
