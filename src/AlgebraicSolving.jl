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
include("algorithms/normal-forms.jl")
include("algorithms/solvers.jl")
#= siggb =#
include("siggb/siggb.jl")
#= examples =#
include("examples/katsura.jl")
include("examples/cyclic.jl")

end # module AlgebraicSolving
