using AlgebraicSolving
using Test

@testset verbose = true "AlgebraicSolving Tests" begin
include("interfaces/nemo.jl")
include("algorithms/groebner-bases.jl")
include("algorithms/normal-forms.jl")
include("algorithms/solvers.jl")
include("algorithms/dimension.jl")
include("algorithms/hilbert.jl")
include("algorithms/decomposition.jl")
include("algorithms/param-curves.jl")
include("examples/katsura.jl")
include("interp/thiele.jl") 
include("interp/newton.jl")
include("interp/cuyt_lee.jl")
end
