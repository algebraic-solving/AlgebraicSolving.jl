module AlgebraicSolving

greet() = print("AlgebraicSolving -- a package for algebraic solving based on msolve")

include("imports.jl")

export PolynomialRing, MPolyRing, QQ, ZZ, FiniteField, GF, characteristic

include("interfaces/abstract-algebra.jl")
include("algorithms/groebner-bases.jl")

end # module AlgebraicSolving
