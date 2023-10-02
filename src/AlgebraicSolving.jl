module AlgebraicSolving

greet() = print("AlgebraicSolving -- a package for algebraic solving based on msolve")

#= io =#
include("imports.jl")
include("exports.jl")
#= types =#
include("types.jl")
#= functionality =#
include("interfaces/nemo.jl")
include("algorithms/groebner-bases.jl")
include("algorithms/solvers.jl")
#= f5 =#
include("f5/f5.jl")
#= examples =#
include("examples/katsura.jl")
include("examples/cyclic.jl")

end # module AlgebraicSolving
