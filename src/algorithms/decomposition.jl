using Logging

@doc Markdown.doc"""
    equidimensional_decomposition(I::Ideal{T}, info_level::Int=0) where {T <: MPolyRingElem}

Given a polynomial ideal `I`, return a list of locally closed sets
`dec` s.t.  each locally closed set in `dec` is equidimensional
(i.e. its Zariski closure has irreducible components of only one
dimension) and s.t. the algebraic set defined by `I` equals the union of
the locally closed sets in `dec`.

When the ground field of `I` is the rational numbers, a multi-modular method is used.

# Arguments
- `I::Ideal{T} where T <: MPolyRingElem`: input ideal.
- `info_level::Int=0`: info level printout: off (`0`, default), details about multi-modular computation (`1`, only relevant when `I` is defined over the rational numbers), computational details (`2`)

# Example
```jldoctest
julia> using AlgebraicSolving

julia> R, (x, y, z) = polynomial_ring(GF(65521), ["x", "y", "z"])
(Multivariate polynomial ring in 3 variables over GF(65521), FqMPolyRingElem[x, y, z])

julia> I = Ideal([x*y - x*z, x*z^2 - x*z, x^2*z - x*z])
FqMPolyRingElem[x*y + 65520*x*z, x*z^2 + 65520*x*z, x^2*z + 65520*x*z]

julia> equidimensional_decomposition(I)
3-element Vector{AlgebraicSolving.LocallyClosedSet{FqMPolyRingElem}}:
 V(x*y + 65520*x*z) \ V(y + 65520*z)
 V(y + 65520*z, x*z^2 + 65520*x*z) \ V(x * (z + 65520))
 V(y + 65520*z, x*z^2 + 65520*x*z, x^2*z + 65520*x*z, z + 65520) \ V(x * z)
```
"""
function equidimensional_decomposition(I::Ideal{T};
                                       info_level::Int=0) where {T <: MPolyRingElem}

    log_level = if info_level == 0
        Warn
    elseif info_level == 1
        INFOONE
    else
        Info
    end
        
    logger = ConsoleLogger(stdout, log_level)
    return _equidimensional_decomposition(I, logger)
end

function _equidimensional_decomposition(I::Ideal{T},
                                        logger::ConsoleLogger) where {T <: FqMPolyRingElem}

    F = I.gens
    Fhom = homogenize(F)
    sort!(Fhom, by = p -> total_degree(p))
    r = ModularRegistry(T[])
    sys_mons, sys_coeffs, basis_ht, char, shift = input_setup(Fhom)
    cells = with_logger(logger) do
        _sig_decomp(sys_mons, sys_coeffs, basis_ht, char, shift, parent(first(Fhom)), r)
    end
    res = LocallyClosedSet{T}[]
    R = parent(I)
    for cell in cells
        append!(res, get_output_cells(cell, R, r))
    end
    return res
end

function _equidimensional_decomposition(I::Ideal{T},
                                        logger::ConsoleLogger) where {T <: QQMPolyRingElem}

    F = I.gens
    Fhom = homogenize(F)
    sort!(Fhom, by = p -> total_degree(p))
    Rhom = parent(first(Fhom))
    r = ReconstructRegistry(Rhom, ReconstructPol[], 1,
                            Int32[], Int32(0))
    cells = LocClosedSet{FqMPolyRingElem}[]
    cnt = 0

    sys_mons = Vector{MonIdx}[]
    sys_coeffs = Vector{Coeff}[]
    basis_ht = initialize_basis_hash_table(Val(ngens(Rhom)))
    char = zero(Coeff)
    shift = zero(Cbuf)
    
    with_logger(logger) do
        while !is_finished(r)
            cnt += 1
            p = Int32(rand_bits_prime(ZZ, 31))

            reset_rand_coeffs!(r.rand_coeffs)
            reset_tracers!(r.tracers)
            new_prime!(r, p)
            S, _ = polynomial_ring(GF(p), ngens(Rhom),
                                   internal_ordering = :degrevlex)

            Fhomp = [reduce_mod_p(f, S) for f in Fhom]
            if cnt > 1
                char = Coeff(p)
                shift = maxshift(char)
                for (i, f) in enumerate(Fhomp)
                    for (j, c) in enumerate(coefficients(f))
                        sys_coeffs[i][j] = Coeff(lift(ZZ, c).d)
                    end
                end
            else
                sys_mons, sys_coeffs, basis_ht, char, shift = input_setup(Fhomp)
            end
                
            cells = _sig_decomp(sys_mons, sys_coeffs, basis_ht, char, shift, parent(first(Fhomp)), r)
            mark_recorded!(r.tracers)
            if ispow2(cnt)
                @logmsg INFOONE "decomposition $cnt with prime $p, $(length(findall(p -> all(p.is_stable), r.pols))) / $(length(r.pols)) polynomials reconstructed"
            end
            isempty(r.pols) && break # catch the case where no splitting happens
        end
        @logmsg INFOONE "$cnt primes used"
        res = LocallyClosedSet{T}[]
        R = parent(I)
        for cell in cells
            append!(res, get_output_cells(cell, R, Fhom, r))
        end
        return res
    end
end

