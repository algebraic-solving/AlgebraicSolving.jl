module Interpolation

using Markdown
using Nemo
import Nemo.Generic: FracFieldElem
using ..Progress
import ..AlgebraicSolving: _map_exponent_vectors, _change_ring

export thiele, newton, cuyt_lee

# Interpolation algorithms
include("thiele.jl")
include("newton.jl")
include("cuyt_lee.jl")

# Applications
include("resultant.jl")

end # module Interpolation
