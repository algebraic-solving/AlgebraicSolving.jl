using AlgebraicSolving
using Test

@testset verbose = true "AlgebraicSolving Tests" begin
include("interfaces/nemo.jl")
include("algorithms/groebner-bases.jl")
include("algorithms/normal-forms.jl")
include("algorithms/solvers.jl")
include("algorithms/dimension.jl")
include("examples/katsura.jl")
end
