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

function f5(sys::Vector{T}; infolevel = Logging.Warn) where {T <: MPolyElem}
    R = first(sys).parent
    Rchar = characteristic(R)

    # check if input is ok
    if Rchar > 2^31
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

    # convert to and initialize our data structures
    nv = nvars(R)
    basis_ht = initialize_basis_hash_table(Val(nv))

    # initialize basis
    sigs = Vector{Monomial{nv}}(undef, init_basis_size)
    sigmasks = Vector{MaskSig}(undef, init_basis_size)
    sigratios = Vector{Monomial{nv}}(undef, init_basis_size)
    lm_masks = Vector{DivMask}(undef, init_basis_size)
    monomials = Vector{Vector{MonIdx}}(undef, init_basis_size)
    coeffs = Vector{Vector{Coeff}}(undef, init_basis_size)
    is_red = Vector{Bool}(undef, init_basis_size)
    syz_sigs = Vector{Monomial{nv}}(undef, init_syz_size)
    syz_masks = Vector{MaskSig}(undef, init_syz_size)
    basis = Basis(sigs, sigmasks, sigratios, lm_masks,
                  monomials, coeffs, is_red,
                  syz_sigs, syz_masks, degs, sysl,
                  init_basis_size, sysl + 1, 0, init_syz_size)

    # initialize pairset
    pairset = LoadVector(Vector{SPair{nv}}(undef, init_pair_size),
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
        @inbounds for j in 1:lf
            m = monomial(SVector{nv}((Exp).(exps[j])))
            eidx = insert_in_hash_table!(basis_ht, m)
            cf = Coeff(cfs[j].data)
            mons[j] = eidx
            coeffs[j] = cf
        end

        # signatures
        sig = (SigIndex(i), one_mon)
        lm_exps = SVector{nv}((Exp).(exps[1]))
        sigr = monomial(-lm_exps)

        # store stuff in basis
        basis.sigs[i] = one_mon
        basis.sigratios[i] = sigr
        basis.monomials[i] = mons
        basis.coefficients[i] = coeffs
        basis.is_red[i] = false

        # add unitvector as pair
        pairset.elems[i] = SPair{nv}(sig, zero_sig, zero(DivMask),
                                     zero(DivMask), i, 0, degs[i])
    end

    # compute divmasks
    dm_one_mon = divmask(one_mon, basis_ht.divmap, basis_ht.ndivbits)
    fill_divmask!(basis_ht)
    for i in 1:sysl
        basis.sigmasks[i] = (SigIndex(i), dm_one_mon)
        pairset.elems[i].top_sig_mask = basis.sigmasks[i][2]
        basis.lm_masks[i] = basis_ht.hashdata[basis.monomials[i][1]].divmask
    end

    # constants for fast arithmetic
    char = Val(Coeff(Rchar.d))
    shift = Val(maxshift(char))

    logger = SimpleLogger(stdout, infolevel)
    with_logger(logger) do
        f5!(basis, pairset, basis_ht, char, shift)
    end

end

function f5!(basis::Basis{N},
             pairset::Pairset,
             basis_ht::MonomialHashtable,
             char::Val{Char},
             shift::Val{Shift}) where {N, Char, Shift}

    while !iszero(pairset.load)
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
    ordr = Base.Sort.ord(isless, identity, false, Base.Sort.Forward)
    sort!(pairset.elems, from, from+sz, def_sort_alg, ordr) 
end
