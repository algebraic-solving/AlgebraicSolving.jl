@doc Markdown.doc"""
    equidimensional_decomposition(I::Ideal{T}, info_level::Int=0) where {T <: MPolyRingElem}

Given a polynomial ideal `I`, return a list of ideals `dec` s.t.
each ideal in `dec` is equidimensional (i.e. has minimal primes
only of one fixed dimension) and s.t. the radical of `I` equals
the intersection of the radicals of the ideals in `dec`.

**Note**: At the moment only ground fields of characteristic `p`, `p` prime, `p < 2^{31}` are supported.

# Arguments
- `I::Ideal{T} where T <: MpolyElem`: input ideal.
- `info_level::Int=0`: info level printout: off (`0`, default), computational details (`1`)

# Example
```jldoctest
julia> using AlgebraicSolving

julia> R, (x, y, z) = polynomial_ring(GF(65521), ["x", "y", "z"])
(Multivariate polynomial ring in 3 variables over GF(65521), FqMPolyRingElem[x, y, z])

julia> I = Ideal([x*y - x*z, x*z^2 - x*z, x^2*z - x*z])
FqMPolyRingElem[x*y + 65520*x*z, x*z^2 + 65520*x*z, x^2*z + 65520*x*z]

julia> equidimensional_decomposition(I)
3-element Vector{Ideal{FqMPolyRingElem}}:
 FqMPolyRingElem[x]
 FqMPolyRingElem[z, y]
 FqMPolyRingElem[z + 65520, y + 65520, x + 65520]
```
"""
function equidimensional_decomposition(I::Ideal{T};
                                       info_level::Int=0) where {T <: MPolyRingElem}

    F = I.gens
    Fhom = homogenize(F)
    sort!(Fhom, by = p -> total_degree(p))
    cells = _sig_decomp(Fhom, info_level = info_level)
    res = typeof(I)[]
    R = parent(I)
    for cell in cells
        for gb in cell.gbs
            gb_dehom = _dehomogenize(gb, R)
            idl = Ideal(gb_dehom)
            idl.gb[0] = gb_dehom
            push!(res, idl)
        end
    end
    return res
end
