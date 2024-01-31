function input_setup(sys::Vector{<:MPolyRingElem})

    if isempty(sys)
        error("Input system is empty.")
    end
    
    R = first(sys).parent
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

    sys_mons = Vector{Vector{MonIdx}}(undef, sysl)
    sys_coeffs = Vector{Vector{Coeff}}(undef, sysl)

    # store initial pols in basis and pairset
    @inbounds for i in 1:sysl
        f = sys[i]
        lf = length(f)

        # gather up monomials and coeffs
        coeffs, mons = convert_to_ht(f, basis_ht, char,
                                     by = eidx -> basis_ht.exponents[eidx],
                                     lt = lt_drl, rev = true)
        sys_mons[i] = copy(mons)
        sys_coeffs[i] = copy(coeffs)
    end

    return sys_mons, sys_coeffs, basis_ht, char, shift
end

function fill_basis!(sys_mons::Vector{Vector{MonIdx}},
                     sys_coeffs::Vector{Vector{Coeff}},
                     basis_ht::MonomialHashtable{N}) where N

    # initialize basis
    sysl = length(sys_mons)
    basis = new_basis(init_basis_size, init_syz_size, sysl, Val(N))

    @inbounds for i in 1:sysl
        s_ind = SigIndex(i)
        mons = sys_mons[i]
        coeffs = sys_coeffs[i]
        lm = basis_ht.exponents[first(mons)]
        lm_mask = divmask(lm, basis_ht.divmap, basis_ht.ndivbits)

        add_input_element!(basis, s_ind,
                           sys_mons[i], sys_coeffs[i],
                           lm_mask, lm)
    end

    return basis
end

function convert_to_pol(R::MPolyRing,
                        exps::Vector{<:Monomial},
                        coeffs::Vector{Coeff})

    ctx = MPolyBuildCtx(R)
    for (e, c) in zip(exps, coeffs)
        push_term!(ctx, base_ring(R)(c), Vector{Int}(e.exps))
    end
    return finish(ctx)
end

function convert_to_ht(f::MPolyRingElem,
                       ht::MonomialHashtable{N},
                       char::Val{Char};
                       normalise::Bool=true,
                       kwargs...) where {N, Char}

    lf = length(f)

    # gather up monomials and coeffs
    exps = collect(exponent_vectors(f))
    cfs = collect(coefficients(f))
    mons = Vector{MonIdx}(undef, lf)
    coeffs = Vector{Coeff}(undef, lf)
    inver = one(Coeff)
    check_enlarge_hashtable!(ht, lf)
    @inbounds for j in 1:lf
        m = monomial(SVector{N}((Exp).(exps[j])))
        eidx = insert_in_hash_table!(ht, m)
        cf = Coeff(lift(ZZ, cfs[j]).d)
        mons[j] = eidx
        coeffs[j] = cf
    end
    s = sortperm(mons; kwargs...)
    mons = mons[s]
    coeffs = coeffs[s]
    if normalise
        inver = inv(coeffs[1], char)
        coeffs[1] = one(Coeff)
        for i in 2:lf
            coeffs[i] = mul(inver, coeffs[i], char)
        end
    end
    return coeffs, mons
end

function new_basis(basis_size, syz_size,
                   input_length, ::Val{N}) where N

    sigs = Vector{Sig{N}}(undef, basis_size)
    sigmasks = Vector{MaskSig}(undef, basis_size)
    sigratios = Vector{Monomial{N}}(undef, basis_size)
    rewrite_nodes = Vector{Vector{Int}}(undef, basis_size+1)
    lm_masks = Vector{DivMask}(undef, basis_size)
    monomials = Vector{Vector{MonIdx}}(undef, basis_size)
    coeffs = Vector{Vector{Coeff}}(undef, basis_size)
    is_red = Vector{Bool}(undef, basis_size)
    mod_rep_known = Vector{Vector{Bool}}(undef, basis_size)
    mod_reps = Vector{Vector{Polynomial}}(undef, basis_size)
    syz_sigs = Vector{Monomial{N}}(undef, syz_size)
    syz_masks = Vector{MaskSig}(undef, syz_size)
    basis = Basis(sigs, sigmasks, sigratios, rewrite_nodes,
                  lm_masks, monomials, coeffs, is_red,
                  mod_rep_known, mod_reps,
                  syz_sigs, syz_masks, Exp[],
                  input_length,
                  init_basis_size, 0, input_length,
                  input_length + 1, 0,
                  init_syz_size)

    # root node
    basis.rewrite_nodes[1] = [-1, -1]

    return basis
end

function init_pairset(::Val{N}) where N
    ps = Pairset(Vector{SPair{N}}(undef, init_pair_size),
                 0, init_pair_size)
    return ps
end

function add_input_element!(basis::Basis{N},
                            ind::SigIndex,
                            mons::Vector{MonIdx},
                            coeffs::Vector{Coeff},
                            lm_divm::DivMask,
                            lm::Monomial) where N

    @inbounds begin
        one_mon = one_monomial(Monomial{N})

        # signature
        sig = (ind, one_mon)

        l = basis.input_load + 1

        # store stuff in basis
        basis.sigs[l] = sig
        basis.sigmasks[l] = (ind, zero(DivMask))
        basis.sigratios[l] = lm
        basis.rewrite_nodes[l+1] = [-1, 1]
        basis.monomials[l] = mons
        basis.coefficients[l] = coeffs
        basis.is_red[l] = false
        push!(basis.degs, lm.deg)
        basis.lm_masks[l] = lm_divm
        basis.input_load += 1

        # add child to rewrite root
        push!(basis.rewrite_nodes[1], l+1)
        basis.rewrite_nodes[1][1] += 1
    end
end

function add_unit_pair!(basis::Basis{N},
                        pairset::Pairset{N},
                        ind::Integer,
                        deg::Exp) where N

    @inbounds basis.is_red[ind] && return
    one_mon = one_monomial(Monomial{N})
    zero_sig = (zero(SigIndex), one_mon)
    @inbounds sig = basis.sigs[ind]
    @inbounds msk = basis.sigmasks[ind]
    pairset.elems[pairset.load+1] = SPair{N}(sig, zero_sig, mask(msk),
                                             zero(DivMask), Int(ind),
                                             0, deg)
    pairset.load += 1
    return
end
