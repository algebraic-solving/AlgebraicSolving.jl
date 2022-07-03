using AlgebraicSolving
using Test

@testset verbose = true "AlgebraicSolving Tests" begin
include("interfaces/abstract-algebra.jl")
include("algorithms/groebner-bases.jl")
include("algorithms/solvers.jl")
include("examples/katsura.jl")
end
