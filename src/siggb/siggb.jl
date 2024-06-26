# sizes for initialization
const init_ht_size = 2^17
const init_basis_size = 10000
const init_syz_size = 1000
const init_pair_size = 10000
# default sorting alg
const def_sort_alg = Base.Sort.DEFAULT_UNSTABLE
include("typedefs.jl")
include("monomials.jl")
include("hashtable.jl")
include("update.jl")
include("symbolic_pp.jl")
include("linear_algebra.jl")

@doc Markdown.doc"""
    sig_groebner_basis(sys::Vector{T}; info_level::Int = 0, degbound::Int = 0) where {T <: MPolyRingElem}

Compute a Signature Gröbner basis of the sequence `sys` w.r.t. to the
degree reverse lexicographical monomial ordering and the degree
position-over-term ordering induced by `sys`. The output is a vector
of `Tuple{Tuple{Int64, T}, T}` where the first element indicates the
signature and the second the underlying polynomial.

**Note**: At the moment only ground fields of characteristic `p`, `p` prime, `p < 2^{31}` are supported.
**Note**: The input generators must be homogeneous.
**Note**: The algorithms behaviour may depend heavily on how the elements in `sys` are sorted.

# Arguments
- `sys::Vector{T} where T <: MpolyElem`: input generators.
- `info_level::Int=0`: info level printout: off (`0`, default), computational details (`1`)
- `degbound::Int=0`: degree bound for Gröbner basis computation, compute a full Gröbner basis if `0` (default) or only up to degree `d`.

# Example
```jldoctest
julia> using AlgebraicSolving

julia> R, vars = polynomial_ring(GF(17), ["x$i" for i in 1:4])
(Multivariate polynomial ring in 4 variables over GF(17), FqMPolyRingElem[x1, x2, x3, x4])

julia> F = AlgebraicSolving.cyclic(R)
FqMPolyRingElem[x1 + x2 + x3 + x4, x1*x2 + x1*x4 + x2*x3 + x3*x4, x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4, x1*x2*x3*x4 + 16]

julia> Fhom = AlgebraicSolving._homogenize(F.gens)
4-element Vector{FqMPolyRingElem}:
 x1 + x2 + x3 + x4
 x1*x2 + x2*x3 + x1*x4 + x3*x4
 x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4
 x1*x2*x3*x4 + 16*x5^4

julia> sig_groebner_basis(Fhom)
7-element Vector{Tuple{Tuple{Int64, FqMPolyRingElem}, FqMPolyRingElem}}:
 ((1, 1), x1 + x2 + x3 + x4)
 ((2, 1), x2^2 + 2*x2*x4 + x4^2)
 ((3, 1), x2*x3^2 + x3^2*x4 + 16*x2*x4^2 + 16*x4^3)
 ((4, 1), x2*x3*x4^2 + x3^2*x4^2 + 16*x2*x4^3 + x3*x4^3 + 16*x4^4 + 16*x5^4)
 ((4, x3), x3^3*x4^2 + x3^2*x4^3 + 16*x3*x5^4 + 16*x4*x5^4)
 ((4, x2), x2*x4^4 + x4^5 + 16*x2*x5^4 + 16*x4*x5^4)
 ((4, x2*x3), x3^2*x4^4 + x2*x3*x5^4 + 16*x2*x4*x5^4 + x3*x4*x5^4 + 15*x4^2*x5^4)
```
"""
function sig_groebner_basis(sys::Vector{T}; info_level::Int=0, degbound::Int=0) where {T <: MPolyRingElem}
    R = parent(first(sys))
    Rchar = characteristic(R)

    # check if input is ok
    if Rchar > 2^31 || iszero(Rchar)
        error("At the moment we only support finite fields up to prime characteristic < 2^31.")
    end
    sysl = length(sys)
    degs = Vector{Exp}(undef, sysl)
    @inbounds for (i, f) in enumerate(sys)
        deg = total_degree(f)
        if deg > typemax(Exp)
            error("input degrees too large.")
        end
        degs[i] = Exp(deg)
        for m in exponent_vectors(f)
            if sum(m) != deg
                error("input system must be homogeneous.")
            end
        end
    end

    # constants for fast arithmetic
    char = Val(Coeff(Rchar.d))
    shift = Val(maxshift(char))

    # convert to and initialize our data structures
    nv = nvars(R)
    basis_ht = initialize_basis_hash_table(Val(nv))

    # initialize basis
    sigs = Vector{Sig{nv}}(undef, init_basis_size)
    sigmasks = Vector{MaskSig}(undef, init_basis_size)
    sigratios = Vector{Monomial{nv}}(undef, init_basis_size)
    rewrite_nodes = Vector{Vector{Int}}(undef, init_basis_size+1)
    lm_masks = Vector{DivMask}(undef, init_basis_size)
    monomials = Vector{Vector{MonIdx}}(undef, init_basis_size)
    coeffs = Vector{Vector{Coeff}}(undef, init_basis_size)
    is_red = Vector{Bool}(undef, init_basis_size)
    syz_sigs = Vector{Monomial{nv}}(undef, init_syz_size)
    syz_masks = Vector{MaskSig}(undef, init_syz_size)
    basis = Basis(sigs, sigmasks, sigratios, rewrite_nodes,
                  lm_masks, monomials, coeffs, is_red,
                  syz_sigs, syz_masks, degs, sysl,
                  init_basis_size, sysl + 1, 0, init_syz_size)

    # root node
    basis.rewrite_nodes[1] = vcat([length(sys)-1, -1],
                                  collect(2:length(sys)+1))

    # initialize pairset
    pairset = Pairset{nv}(Vector{SPair{nv}}(undef, init_pair_size),
                          sysl,
                          init_pair_size)

    one_mon = monomial(SVector{nv}(zeros(Exp, nv)))
    zero_sig = (zero(SigIndex), one_mon)
    # store initial pols in basis and pairset
    @inbounds for i in 1:sysl
        f = sys[i]
        lf = length(f)

        # gather up monomials and coeffs
        exps = collect(exponent_vectors(f))
        cfs = collect(coefficients(f))
        mons = Vector{MonIdx}(undef, lf)
        coeffs = Vector{Coeff}(undef, lf)
        inver = one(Coeff)
        @inbounds for j in 1:lf
            m = monomial(SVector{nv}((Exp).(exps[j])))
            eidx = insert_in_hash_table!(basis_ht, m)
            if isone(j)
                inver = inv(Coeff(lift(ZZ, cfs[1])), char)
            end
            cf = isone(j) ? one(Coeff) : mul(inver, Coeff(lift(ZZ, cfs[j])), char)
            mons[j] = eidx
            coeffs[j] = cf
        end
        s = sortperm(mons, by = eidx -> basis_ht.exponents[eidx],
                     lt = lt_drl, rev = true)
        @inbounds mons = mons[s]
        @inbounds coeffs = coeffs[s]

        # signatures
        sig = (SigIndex(i), one_mon)
        lm_exps = SVector{nv}((Exp).(exps[1]))
        sigr = monomial(lm_exps)

        # store stuff in basis
        basis.sigs[i] = sig
        basis.sigratios[i] = sigr
        basis.rewrite_nodes[i+1] = [-1, 1]
        basis.monomials[i] = mons
        basis.coefficients[i] = coeffs
        basis.is_red[i] = false

        # add unitvector as pair
        pairset.elems[i] = SPair{nv}(sig, zero_sig, zero(DivMask),
                                     zero(DivMask), i, 0, degs[i])
    end

    # compute divmasks
    fill_divmask!(basis_ht)
    dm_one_mon = divmask(one_mon, basis_ht.divmap, basis_ht.ndivbits)
    @inbounds for i in 1:sysl
        basis.sigmasks[i] = (SigIndex(i), dm_one_mon)
        pairset.elems[i].top_sig_mask = basis.sigmasks[i][2]
        basis.lm_masks[i] = basis_ht.hashdata[basis.monomials[i][1]].divmask
    end

    logger = ConsoleLogger(stdout, info_level == 0 ? Warn : Info)
    with_logger(logger) do
        siggb!(basis, pairset, basis_ht, char, shift, degbound = degbound)
    end

    # output
    eltp = typeof(first(sys))
    outp = Tuple{Tuple{Int, eltp}, eltp}[]
    @inbounds for i in basis.basis_offset:basis.basis_load
        exps = [basis_ht.exponents[m].exps for m in basis.monomials[i]]
        ctx = MPolyBuildCtx(R)
        for (e, c) in zip(exps, basis.coefficients[i])
            push_term!(ctx, coefficient_ring(R)(c), Vector{Int}(e))
        end
        pol = finish(ctx)

        s = basis.sigs[i]
        ctx = MPolyBuildCtx(R)
        push_term!(ctx, one(coefficient_ring(R)), Vector{Int}(monomial(s).exps))
        sig = (Int(index(s)), finish(ctx))

        push!(outp, (sig, pol))
    end
    return outp
end

function siggb!(basis::Basis{N},
                pairset::Pairset,
                basis_ht::MonomialHashtable,
                char::Val{Char},
                shift::Val{Shift};
                degbound = 0) where {N, Char, Shift}

    while !iszero(pairset.load)
        if !iszero(degbound) && first(pairset.elems).deg > degbound
            break
        end
	matrix = initialize_matrix(Val(N))
        symbol_ht = initialize_secondary_hash_table(basis_ht)

        select_normal!(pairset, basis, matrix, basis_ht, symbol_ht)
        symbolic_pp!(basis, matrix, basis_ht, symbol_ht)
        finalize_matrix!(matrix, symbol_ht)
        echelonize!(matrix, char, shift)

        update_basis!(basis, matrix, pairset, symbol_ht, basis_ht)
    end
end


# miscallaneous helper functions
function sort_pairset_by_degree!(pairset::Pairset, from::Int, sz::Int)
    ordr = Base.Sort.ord(isless, p -> p.deg, false, Base.Sort.Forward)
    sort!(pairset.elems, from, from+sz, def_sort_alg, ordr) 
end

# homogenize w.r.t. the last variable
function _homogenize(F::Vector{P}) where {P <: MPolyRingElem}
    R = parent(first(F))
    S, vars = polynomial_ring(base_ring(R), ["x$i" for i in 1:nvars(R)+1],
                              internal_ordering = :degrevlex)
    res = typeof(first(F))[]
    for f in F
        ctx = MPolyBuildCtx(S)
        d = total_degree(f)
        for (e, c) in zip(exponent_vectors(f), coefficients(f))
            enew = push!(e, d - sum(e))
            push_term!(ctx, c, e)
        end
        push!(res, finish(ctx))
    end
    return res
end

# test against msolve
function _is_gb(gb::Vector{Tuple{Tuple{Int, P}, P}}) where {P <: MPolyRingElem}
    gb_pols = [p[2] for p in gb]
    gb_msolve = groebner_basis(Ideal(gb_pols), complete_reduction = true)
    
    lms_gb = (Nemo.leading_monomial).(gb_pols)
    lms_msolve = (Nemo.leading_monomial).(gb_msolve)
    res1 = all(u -> any(v -> divides(u, v)[1], lms_gb), lms_msolve)
    res2 = all(u -> any(v -> divides(u, v)[1], lms_msolve), lms_gb)
    return res1 && res2
end
