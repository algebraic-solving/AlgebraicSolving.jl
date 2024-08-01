import msolve_jll: libneogb, libmsolve

export groebner_basis, eliminate

@doc Markdown.doc"""
    eliminate(I::Ideal{T} where T <: MPolyRingElem, eliminate::Int,  <keyword arguments>)

Compute a Groebner basis of the ideal `I` w.r.t. to the product monomial ordering defined by two blocks
w.r.t. the degree reverse lexicographical monomial ordering using Faugère's F4 algorithm. Hereby the first block includes
the first `eliminate` variables.

At the moment the underlying algorithm is based on variants of Faugère's F4 Algorithm.

**Note**: At the moment only ground fields of characteristic `p`, `p` prime, `p < 2^{31}` and the rationals are supported.

# Arguments
- `I::Ideal{T} where T <: MPolyRingElem`: input generators.
- `eliminate::Int=0`: size of first block of variables to be eliminated.
- `intersect::Bool=true`: compute the `eliminate`-th elimination ideal.
- `initial_hts::Int=17`: initial hash table size `log_2`.
- `nr_thrds::Int=1`: number of threads for parallel linear algebra.
- `max_nr_pairs::Int=0`: maximal number of pairs per matrix, only bounded by minimal degree if `0`.
- `la_option::Int=2`: linear algebra option: exact sparse-dense (`1`), exact sparse (`2`, default), probabilistic sparse-dense (`42`), probabilistic sparse(`44`).
- `complete_reduction::Bool=true`: compute a reduced Gröbner basis for `I`.
- `normalize::Bool=false`: normalize generators of Gröbner basis for `I`, only applicable when working over the rationals.
- `truncate_lifting::Int=0`: truncates the lifting process to given number of elements, only applicable when working over the rationals.
- `info_level::Int=0`: info level printout: off (`0`, default), summary (`1`), detailed (`2`).

# Examples
```jldoctest
julia> using AlgebraicSolving

julia> R, (x,y,z) = polynomial_ring(GF(101),["x","y","z"], internal_ordering=:degrevlex)
(Multivariate polynomial ring in 3 variables over GF(101), FqMPolyRingElem[x, y, z])

julia> I = Ideal([x+2*y+2*z-1, x^2+2*y^2+2*z^2-x, 2*x*y+2*y*z-y])
FqMPolyRingElem[x + 2*y + 2*z + 100, x^2 + 2*y^2 + 2*z^2 + 100*x, 2*x*y + 2*y*z + 100*y]

julia> eliminate(I, 2)
1-element Vector{FqMPolyRingElem}:
 z^4 + 38*z^3 + 95*z^2 + 95*z
```
"""
function eliminate(
        I::Ideal{T} where T <: MPolyRingElem,
        eliminate::Int;
        intersect::Bool=true,
        initial_hts::Int=17,
        nr_thrds::Int=1,
        max_nr_pairs::Int=0,
        la_option::Int=2,
        complete_reduction::Bool=true,
        normalize::Bool=false,
        truncate_lifting::Int=0,
        info_level::Int=0
        )
    if eliminate <= 0
        error("Number of variables to be eliminated is <= 0.")
    else
        return groebner_basis(I, initial_hts=initial_hts,nr_thrds=nr_thrds,
                              max_nr_pairs=max_nr_pairs, la_option=la_option,
                              eliminate=eliminate, intersect=intersect,
                              complete_reduction=complete_reduction,
                              normalize = normalize,
                              truncate_lifting = truncate_lifting,
                              info_level=info_level)
    end
end

@doc Markdown.doc"""
    groebner_basis(I::Ideal{T} where T <: MPolyRingElem, <keyword arguments>)

Compute a Groebner basis of the ideal `I` w.r.t. to the degree reverse lexicographical monomial ordering using Faugère's F4 algorithm.
At the moment the underlying algorithm is based on variants of Faugère's F4 Algorithm.

**Note**: At the moment only ground fields of characteristic `p`, `p` prime, `p < 2^{31}` and the rationals are supported.

# Arguments
- `I::Ideal{T} where T <: MPolyRingElem`: input generators.
- `initial_hts::Int=17`: initial hash table size `log_2`.
- `nr_thrds::Int=1`: number of threads for parallel linear algebra.
- `max_nr_pairs::Int=0`: maximal number of pairs per matrix, only bounded by minimal degree if `0`.
- `la_option::Int=2`: linear algebra option: exact sparse-dense (`1`), exact sparse (`2`, default), probabilistic sparse-dense (`42`), probabilistic sparse(`44`).
- `eliminate::Int=0`: size of first block of variables to be eliminated.
- `intersect::Bool=true`: compute the `eliminate`-th elimination ideal.
- `complete_reduction::Bool=true`: compute a reduced Gröbner basis for `I`.
- `normalize::Bool=false`: normalize generators of Gröbner basis for `I`, only applicable when working over the rationals.
- `truncate_lifting::Int=0`: truncates the lifting process to given number of elements, only applicable when working over the rationals.
- `info_level::Int=0`: info level printout: off (`0`, default), summary (`1`), detailed (`2`).

# Examples
```jldoctest
julia> using AlgebraicSolving

julia> R, (x,y,z) = polynomial_ring(GF(101),["x","y","z"], internal_ordering=:degrevlex)
(Multivariate polynomial ring in 3 variables over GF(101), FqMPolyRingElem[x, y, z])

julia> I = Ideal([x+2*y+2*z-1, x^2+2*y^2+2*z^2-x, 2*x*y+2*y*z-y])
FqMPolyRingElem[x + 2*y + 2*z + 100, x^2 + 2*y^2 + 2*z^2 + 100*x, 2*x*y + 2*y*z + 100*y]

julia> groebner_basis(I)
4-element Vector{FqMPolyRingElem}:
 x + 2*y + 2*z + 100
 y*z + 82*z^2 + 10*y + 40*z
 y^2 + 60*z^2 + 20*y + 81*z
 z^3 + 28*z^2 + 64*y + 13*z

julia> groebner_basis(I, eliminate=2)
1-element Vector{FqMPolyRingElem}:
 z^4 + 38*z^3 + 95*z^2 + 95*z
```
"""
function groebner_basis(
        I::Ideal{T} where T <: MPolyRingElem;
        initial_hts::Int=17,
        nr_thrds::Int=1,
        max_nr_pairs::Int=0,
        la_option::Int=2,
        eliminate::Int=0,
        intersect::Bool=true,
        complete_reduction::Bool=true,
        normalize::Bool=false,
        truncate_lifting::Int=0,
        info_level::Int=0
        )

    return get!(I.gb, eliminate) do
        _core_groebner_basis(I, initial_hts = initial_hts, nr_thrds = nr_thrds, 
                             max_nr_pairs = max_nr_pairs, la_option = la_option,
                             eliminate = eliminate, intersect = intersect,
                             complete_reduction = complete_reduction,
                             normalize = normalize,
                             truncate_lifting = truncate_lifting,
                             info_level = info_level)
    end
end

function _core_groebner_basis(
        I::Ideal{T} where T <: MPolyRingElem;
        initial_hts::Int=17,
        nr_thrds::Int=1,
        max_nr_pairs::Int=0,
        la_option::Int=2,
        eliminate::Int=0,
        intersect::Bool=true,
        complete_reduction::Bool=true,
        normalize::Bool=false,
        truncate_lifting::Int=0,
        info_level::Int=0
        )

    F = I.gens
    R = first(F).parent

    if F == repeat([R(0)], length(F))
        I.gb[eliminate] = F
        return F
    end
    nr_vars     = nvars(R)
    nr_gens     = length(F)
    field_char  = Int(characteristic(R))

    mon_order       = 0
    elim_block_size = eliminate
    if elim_block_size >= nr_vars
        error("Number of variables to be eliminated is bigger than number of variables in ring.")
    end
    reduce_gb       = Int(complete_reduction)

    # convert ideal to flattened arrays of ints

    if !(field_char == 0)
        if !(is_probable_prime(field_char))
            error("At the moment we only support finite fields or the rationals.")
        end
    end

    # nr_gens might change if F contains zero polynomials
    lens, cfs, exps, nr_gens = _convert_to_msolve(F)

    gb_ld  = Ref(Cint(0))
    gb_len = Ref(Ptr{Cint}(0))
    gb_exp = Ref(Ptr{Cint}(0))
    gb_cf  = Ref(Ptr{Cvoid}(0))

    if field_char == 0
        nr_terms  = ccall((:export_groebner_qq, libmsolve), Int,
            (Ptr{Nothing}, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cvoid}},
            Ptr{Cint}, Ptr{Cint}, Ptr{Cvoid}, Cint, Cint, Cint, Cint, Cint, Cint,
            Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint),
            cglobal(:jl_malloc), gb_ld, gb_len, gb_exp, gb_cf, lens, exps, cfs,
            field_char, mon_order, elim_block_size, nr_vars, nr_gens, initial_hts,
            nr_thrds, max_nr_pairs, 0, la_option, reduce_gb, 0, truncate_lifting, info_level)

    else
        nr_terms  = ccall((:export_f4, libneogb), Int,
            (Ptr{Nothing}, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cvoid}},
            Ptr{Cint}, Ptr{Cint}, Ptr{Cvoid}, Cint, Cint, Cint, Cint, Cint, Cint,
            Cint, Cint, Cint, Cint, Cint, Cint, Cint),
            cglobal(:jl_malloc), gb_ld, gb_len, gb_exp, gb_cf, lens, exps, cfs,
            field_char, mon_order, elim_block_size, nr_vars, nr_gens, initial_hts,
            nr_thrds, max_nr_pairs, 0, la_option, reduce_gb, 0, info_level)
    end

    # convert to julia array, also give memory management to julia
    jl_ld   = gb_ld[]
    jl_len  = Base.unsafe_wrap(Array, gb_len[], jl_ld)
    jl_exp  = Base.unsafe_wrap(Array, gb_exp[], nr_terms*nr_vars)

    #  coefficient handling depending on field characteristic
    if field_char == 0
        ptr     = reinterpret(Ptr{BigInt}, gb_cf[])
        jl_cf   = [QQFieldElem(unsafe_load(ptr, i)) for i in 1:nr_terms]
    else
        ptr     = reinterpret(Ptr{Int32}, gb_cf[])
        jl_cf   = Base.unsafe_wrap(Array, ptr, nr_terms)
    end

    #  shall eliminated variables be removed?
    if !intersect
        eliminate = 0
    end

    #  convert to basis
    if field_char == 0
        basis = _convert_rational_array_to_abstract_algebra(
            jl_ld, jl_len, jl_cf, jl_exp, R, normalize, eliminate)
    else
        basis = _convert_finite_field_array_to_abstract_algebra(
            jl_ld, jl_len, jl_cf, jl_exp, R, eliminate)
    end
    ccall((:free_f4_julia_result_data, libneogb), Nothing ,
          (Ptr{Nothing}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}},
           Ptr{Ptr{Cvoid}}, Int, Int),
          cglobal(:jl_free), gb_len, gb_exp, gb_cf, jl_ld, field_char)

    I.gb[eliminate] = basis

    return basis
end
