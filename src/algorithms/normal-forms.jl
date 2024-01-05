import msolve_jll: libneogb

export normal_form

@doc Markdown.doc"""
    normal_form(
        f::T,
        I::Ideal{T};
        nr_thrds::Int=1,
        info_level::Int=0
        ) where T <: MPolyRingElem

Compute the normal forms of the elements of `F` w.r.t. a degree reverse
lexicographical Gröbner basis of `I`.

**Note:** If `I` has not already cached a degree reverse lexicographical
Gröbner basis, then this one is first computed.

# Arguments
- `F::Vector{T} where T <: MPolyRingElem`: elements to be reduced.
- `I::Ideal{T} where T <: MPolyRingElem`: ideal data to reduce with.
- `nr_thrds::Int=1`: number of threads for parallel linear algebra.
- `info_level::Int=0`: info level printout: off (`0`, default), summary (`1`), detailed (`2`).

# Examples
```jldoctest
julia> using AlgebraicSolving

julia> R, (x,y) = polynomial_ring(GF(101),["x","y"])
(Multivariate polynomial ring in 2 variables over GF(101), FqMPolyRingElem[x, y])

julia> I = Ideal([y*x^3+12-y, x+y])
FqMPolyRingElem[x^3*y + 100*y + 12, x + y]

julia> f = 2*x^2+7*x*y
2*x^2 + 7*x*y

julia> normal_form(f, I)
96*y^2
```
"""
function normal_form(
        f::T,
        I::Ideal{T};
        nr_thrds::Int=1,
        info_level::Int=0
        ) where T <: MPolyRingElem
    nf = _core_normal_form([f], I; nr_thrds, info_level)
    return nf[1]
end

@doc Markdown.doc"""
    normal_form(
        F::Vector{T},
        I::Ideal{T};
        nr_thrds::Int=1,
        info_level::Int=0
        ) where T <: MPolyRingElem

Compute the normal forms of the elements of `F` w.r.t. a degree reverse
lexicographical Gröbner basis of `I`.

**Note:** If `I` has not already cached a degree reverse lexicographical
Gröbner basis, then this one is first computed.

# Arguments
- `F::Vector{T} where T <: MPolyRingElem`: elements to be reduced.
- `I::Ideal{T} where T <: MPolyRingElem`: ideal data to reduce with.
- `nr_thrds::Int=1`: number of threads for parallel linear algebra.
- `info_level::Int=0`: info level printout: off (`0`, default), summary (`1`), detailed (`2`).

# Examples
```jldoctest
julia> using AlgebraicSolving

julia> R, (x,y) = polynomial_ring(GF(101),["x","y"])
(Multivariate polynomial ring in 2 variables over GF(101), FqMPolyRingElem[x, y])

julia> I = Ideal([y*x^3+12-y, x+y])
FqMPolyRingElem[x^3*y + 100*y + 12, x + y]

julia> F = [2*x^2+7*x*y, x+y]
2-element Vector{FqMPolyRingElem}:
 2*x^2 + 7*x*y
 x + y

julia> normal_form(F,I)
2-element Vector{FqMPolyRingElem}:
 96*y^2
 0
```
"""
function normal_form(
        F::Vector{T},
        I::Ideal{T};
        nr_thrds::Int=1,
        info_level::Int=0
        ) where T <: MPolyRingElem
    return _core_normal_form(F, I; nr_thrds, info_level)
end

function _core_normal_form(
        F::Vector{T},
        I::Ideal{T};
        nr_thrds::Int=1,
        info_level::Int=0
        ) where T <: MPolyRingElem

    
    if (length(F) == 0 || length(I.gens) == 0)
        error("Input data not valid.")
    end

    R = first(F).parent
    nr_vars     = nvars(R)
    field_char  = Int(characteristic(R))

    if !(is_probable_prime(field_char))
        error("At the moment we only supports finite fields.")
    end

    #= first get a degree reverse lexicographical Gröbner basis for I =#
    if !haskey(I.gb, 0)
        if field_char > 2^17
            groebner_basis(I, eliminate = 0, la_option = 44, info_level = info_level)
        else
            groebner_basis(I, eliminate = 0, la_option = 2, info_level = info_level)
        end
    end
    G = I.gb[0]

    if G == [R(0)]
        return F
    end

    tbr_nr_gens = length(F)
    bs_nr_gens  = length(G)
    is_gb       = 1

    # convert ideal to flattened arrays of ints
    tbr_lens, tbr_cfs, tbr_exps = _convert_to_msolve(F)
    bs_lens, bs_cfs, bs_exps    = _convert_to_msolve(G)

    nf_ld  = Ref(Cint(0))
    nf_len = Ref(Ptr{Cint}(0))
    nf_exp = Ref(Ptr{Cint}(0))
    nf_cf  = Ref(Ptr{Cvoid}(0))

    nr_terms  = ccall((:export_nf, libneogb), Int,
        (Ptr{Nothing}, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cvoid}},
        Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cvoid},
        Cint, Cint, Cint, Cint, Cint, Cint, Cint),
        cglobal(:jl_malloc), nf_ld, nf_len, nf_exp, nf_cf, tbr_nr_gens, tbr_lens, tbr_exps,
        tbr_cfs, bs_nr_gens, bs_lens, bs_exps, bs_cfs, field_char, 0, 0, nr_vars, is_gb,
        nr_thrds, info_level)

    # convert to julia array, also give memory management to julia
    jl_ld   = nf_ld[]
    jl_len  = Base.unsafe_wrap(Array, nf_len[], jl_ld)
    jl_exp  = Base.unsafe_wrap(Array, nf_exp[], nr_terms*nr_vars)
    ptr     = reinterpret(Ptr{Int32}, nf_cf[])
    jl_cf   = Base.unsafe_wrap(Array, ptr, nr_terms)

    basis = _convert_finite_field_array_to_abstract_algebra(
                jl_ld, jl_len, jl_cf, jl_exp, R, 0)

    ccall((:free_f4_julia_result_data, libneogb), Nothing ,
          (Ptr{Nothing}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}},
           Ptr{Ptr{Cvoid}}, Int, Int),
          cglobal(:jl_free), nf_len, nf_exp, nf_cf, jl_ld, field_char)

    return basis
end

