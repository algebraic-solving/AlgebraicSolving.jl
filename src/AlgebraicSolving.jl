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
include("algorithms/dimension.jl")
include("algorithms/decomposition.jl")
include("algorithms/curve-param.jl")
include("algorithms/hilbert.jl")
include("algorithms/sampling.jl")
#= param-ideal =#
include("algorithms/param-ideal/groebner-bases.jl")
include("algorithms/param-ideal/multiplication-matrices.jl")
include("algorithms/param-ideal/hermite-matrices.jl")
include("algorithms/param-ideal/parametrizations.jl")
include("algorithms/param-ideal/sign-determination.jl")
include("algorithms/param-ideal/semialgebraic-set.jl")
include("algorithms/param-ideal/real-root-classification.jl")
include("algorithms/param-ideal/solve-lmi.jl")
#= siggb =#
include("siggb/siggb.jl")
#= progress =#
include("progress/main.jl")
#= interp =#
include("interp/main.jl")
#= connectivity =#
include("connectivity/connectcurves.jl")
include("connectivity/roadmap.jl")
include("connectivity/polar.jl")
#= examples =#
include("examples/katsura.jl")
include("examples/cyclic.jl")
include("examples/eco.jl")


end # module AlgebraicSolving
