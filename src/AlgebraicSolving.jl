module AlgebraicSolving

greet() = print("AlgebraicSolving -- a package for algebraic solving based on msolve")

include("imports.jl")

export PolynomialRing, MPolyRing, FiniteField, GF, characteristic, ZZ, QQ

include("interfaces/abstract-algebra.jl")
include("algorithms/groebner-bases.jl")
include("algorithms/solvers.jl")
include("examples/katsura.jl")

end # module AlgebraicSolving
