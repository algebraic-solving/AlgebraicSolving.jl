import msolve_jll: libneogb

export groebner_basis, eliminate

@doc Markdown.doc"""
    eliminate(I::Ideal{T} where T <: MPolyRingElem, eliminate::Int,  <keyword arguments>)

Compute a Groebner basis of the ideal `I` w.r.t. to the product monomial ordering defined by two blocks
w.r.t. the degree reverse lexicographical monomial ordering using Faugère's F4 algorithm. Hereby the first block includes
the first `eliminate` variables.

At the moment the underlying algorithm is based on variants of Faugère's F4 Algorithm.

**Note**: At the moment only ground fields of characteristic `p`, `p` prime, `p < 2^{31}` are supported.

# Arguments
- `I::Ideal{T} where T <: MPolyRingElem`: input generators.
- `initial_hts::Int=17`: initial hash table size `log_2`.
- `nr_thrds::Int=1`: number of threads for parallel linear algebra.
- `max_nr_pairs::Int=0`: maximal number of pairs per matrix, only bounded by minimal degree if `0`.
- `la_option::Int=2`: linear algebra option: exact sparse-dense (`1`), exact sparse (`2`, default), probabilistic sparse-dense (`42`), probabilistic sparse(`44`).
- `complete_reduction::Bool=true`: compute a reduced Gröbner basis for `I`
- `info_level::Int=0`: info level printout: off (`0`, default), summary (`1`), detailed (`2`).

# Examples
```jldoctest
julia> using AlgebraicSolving

julia> R, (x,y,z) = polynomial_ring(GF(101),["x","y","z"], ordering=:degrevlex)
(Multivariate polynomial ring in 3 variables over GF(101), fpMPolyRingElem[x, y, z])

julia> I = Ideal([x+2*y+2*z-1, x^2+2*y^2+2*z^2-x, 2*x*y+2*y*z-y])
fpMPolyRingElem[x + 2*y + 2*z + 100, x^2 + 2*y^2 + 2*z^2 + 100*x, 2*x*y + 2*y*z + 100*y]

julia> eliminate(I, 2)
1-element Vector{fpMPolyRingElem}:
 z^4 + 38*z^3 + 95*z^2 + 95*z
```
"""
function eliminate(
        I::Ideal{T} where T <: MPolyRingElem,
        eliminate::Int;
        initial_hts::Int=17,
        nr_thrds::Int=1,
        max_nr_pairs::Int=0,
        la_option::Int=2,
        complete_reduction::Bool=true,
        info_level::Int=0
        )
    if eliminate <= 0
        error("Number of variables to be eliminated is <= 0.")
    else
        return groebner_basis(I, initial_hts=initial_hts, nr_thrds=nr_thrds,
                              max_nr_pairs=max_nr_pairs, la_option=la_option,
                              eliminate=eliminate,
                              complete_reduction=complete_reduction,
                              info_level=info_level)
    end
end

@doc Markdown.doc"""
    groebner_basis(I::Ideal{T} where T <: MPolyRingElem, <keyword arguments>)

Compute a Groebner basis of the ideal `I` w.r.t. to the degree reverse lexicographical monomial ordering using Faugère's F4 algorithm.
At the moment the underlying algorithm is based on variants of Faugère's F4 Algorithm.

**Note**: At the moment only ground fields of characteristic `p`, `p` prime, `p < 2^{31}` are supported.

# Arguments
- `I::Ideal{T} where T <: MPolyRingElem`: input generators.
- `initial_hts::Int=17`: initial hash table size `log_2`.
- `nr_thrds::Int=1`: number of threads for parallel linear algebra.
- `max_nr_pairs::Int=0`: maximal number of pairs per matrix, only bounded by minimal degree if `0`.
- `la_option::Int=2`: linear algebra option: exact sparse-dense (`1`), exact sparse (`2`, default), probabilistic sparse-dense (`42`), probabilistic sparse(`44`).
- `eliminate::Int=0`: size of first block of variables to be eliminated.
- `complete_reduction::Bool=true`: compute a reduced Gröbner basis for `I`
- `info_level::Int=0`: info level printout: off (`0`, default), summary (`1`), detailed (`2`).

# Examples
```jldoctest
julia> using AlgebraicSolving

julia> R, (x,y,z) = polynomial_ring(GF(101),["x","y","z"], ordering=:degrevlex)
(Multivariate polynomial ring in 3 variables over GF(101), fpMPolyRingElem[x, y, z])

julia> I = Ideal([x+2*y+2*z-1, x^2+2*y^2+2*z^2-x, 2*x*y+2*y*z-y])
fpMPolyRingElem[x + 2*y + 2*z + 100, x^2 + 2*y^2 + 2*z^2 + 100*x, 2*x*y + 2*y*z + 100*y]

julia> groebner_basis(I)
4-element Vector{fpMPolyRingElem}:
 x + 2*y + 2*z + 100
 y*z + 82*z^2 + 10*y + 40*z
 y^2 + 60*z^2 + 20*y + 81*z
 z^3 + 28*z^2 + 64*y + 13*z

julia> groebner_basis(I, eliminate=2)
1-element Vector{fpMPolyRingElem}:
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
        complete_reduction::Bool=true,
        info_level::Int=0
        )

    return get!(I.gb, eliminate) do
        _core_groebner_basis(I, initial_hts = initial_hts, nr_thrds = nr_thrds, 
                             max_nr_pairs = max_nr_pairs, la_option = la_option,
                             eliminate = eliminate,
                             complete_reduction = complete_reduction,
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
        complete_reduction::Bool=true,
        info_level::Int=0
        )

    F = I.gens
    R = first(F).parent
    if F == [R(0)]
        return F
    end
    nr_vars     = nvars(R)
    nr_gens     = length(F)
    field_char  = Int(characteristic(R))

    println(field_char)
    println(I)

    mon_order       = 0
    elim_block_size = eliminate
    reduce_gb       = Int(complete_reduction)

    # convert ideal to flattened arrays of ints
    if !(is_probable_prime(field_char))
        error("At the moment we only supports finite fields.")
    end

    lens, cfs, exps = _convert_to_msolve(F)

    gb_ld  = Ref(Cint(0))
    gb_len = Ref(Ptr{Cint}(0))
    gb_exp = Ref(Ptr{Cint}(0))
    gb_cf  = Ref(Ptr{Cvoid}(0))

    nr_terms  = ccall((:export_f4, libneogb), Int,
        (Ptr{Nothing}, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cvoid}},
        Ptr{Cint}, Ptr{Cint}, Ptr{Cvoid}, Cint, Cint, Cint, Cint, Cint, Cint,
        Cint, Cint, Cint, Cint, Cint, Cint, Cint),
        cglobal(:jl_malloc), gb_ld, gb_len, gb_exp, gb_cf, lens, exps, cfs,
        field_char, mon_order, elim_block_size, nr_vars, nr_gens, initial_hts,
        nr_thrds, max_nr_pairs, 0, la_option, reduce_gb, 0, info_level)

    # convert to julia array, also give memory management to julia
    jl_ld   = gb_ld[]
    jl_len  = Base.unsafe_wrap(Array, gb_len[], jl_ld)
    jl_exp  = Base.unsafe_wrap(Array, gb_exp[], nr_terms*nr_vars)
    ptr     = reinterpret(Ptr{Int32}, gb_cf[])
    jl_cf   = Base.unsafe_wrap(Array, ptr, nr_terms)

    basis = _convert_finite_field_array_to_abstract_algebra(
                jl_ld, jl_len, jl_cf, jl_exp, R, eliminate)


    ccall((:free_f4_julia_result_data, libneogb), Nothing ,
          (Ptr{Nothing}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}},
           Ptr{Ptr{Cvoid}}, Int, Int),
          cglobal(:jl_free), gb_len, gb_exp, gb_cf, jl_ld, field_char)

    I.gb[eliminate] = basis

    return basis
end
