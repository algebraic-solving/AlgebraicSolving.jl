function input_setup(sys::Vector{<:MPolyRingElem}, mod_ord::Symbol=:DPOT)

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
        if mod_ord == :DPOT
            degs[i] = Exp(deg)
            for m in exponent_vectors(f)
                if sum(m) != deg
                    error("input system must be homogeneous.")
                end
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

function fill_structs!(sys_mons::Vector{Vector{MonIdx}},
                       sys_coeffs::Vector{Vector{Coeff}},
                       basis_ht::MonomialHashtable{N};
                       sysl::Int=length(sys_mons),
                       def_tg::Symbol=:seq,
                       trace::Val{Bl}=Val(false)) where {N, Bl}

    # initialize basis
    basis = new_basis(init_basis_size, init_syz_size, sysl, Val(N))
    pairset = init_pairset(Val(N))
    tags = Tags()
    ind_order = IndOrder(SigIndex[],
                         Dict{Tuple{SigIndex, SigIndex}, Bool}(),
                         zero(SigIndex))
    tr = Bl ? new_tracer() : NoTracer()

    @inbounds for i in 1:length(sys_mons)
        mons = sys_mons[i]
        coeffs = sys_coeffs[i]
        lm = basis_ht.exponents[first(mons)]
        lm_mask = divmask(lm, basis_ht.divmap, basis_ht.ndivbits)

        add_new_sequence_element!(basis, basis_ht, tr, coeffs, mons,
                                  ind_order, SigIndex(i),
                                  pairset, tags, new_tg = def_tg)
    end

    return basis, pairset, tags, ind_order, tr
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
    sort_poly!((coeffs, mons); kwargs...)
    if normalise
        normalize_cfs!(coeffs, char)
    end
    return coeffs, mons
end

function add_input_element!(basis::Basis{N},
                            mons::Vector{MonIdx},
                            coeffs::Vector{Coeff},
                            lm_divm::DivMask,
                            lm::Monomial) where N

    @inbounds begin
        one_mon = one_monomial(Monomial{N})

        l = basis.input_load + 1

        # signature
        sig = (SigIndex(l), one_mon)

        # store stuff in basis
        basis.sigs[l] = sig
        basis.sigmasks[l] = (SigIndex(l), zero(DivMask))
        basis.sigratios[l] = lm
        basis.rewrite_nodes[l+1] = [-1, 1]
        basis.monomials[l] = mons
        basis.coefficients[l] = coeffs
        basis.is_red[l] = false
        push!(basis.degs, lm.deg)
        basis.lm_masks[l] = lm_divm
        basis.input_load += 1
        for i in basis.basis_offset:basis.basis_load
            push!(basis.mod_rep_known[i], false)
            resize!(basis.mod_reps[i], l)
        end

        # add child to rewrite root
        push!(basis.rewrite_nodes[1], l+1)
        basis.rewrite_nodes[1][1] += 1

        return index(sig)
    end
end

function add_new_sequence_element!(basis::Basis{N},
                                   basis_ht::MonomialHashtable{N},
                                   tr::Tracer,
                                   coeffs::Vector{Coeff},
                                   mons::Vector{MonIdx},
                                   ind_ord::IndOrder,
                                   ord_ind::SigIndex,
                                   pairset::Pairset{N},
                                   tags::Tags;
                                   new_tg::Symbol=:split) where N

    ins_index!(ind_ord, ord_ind)

    # add cofactor to input sequence
    make_room_new_input_el!(basis, pairset, tr)

    lm = basis_ht.exponents[first(mons)]

    lm_divm = divmask(lm, basis_ht.divmap, basis_ht.ndivbits)
    s_ind = add_input_element!(basis, mons,
                               coeffs, lm_divm, lm)

    # update tags
    tags[s_ind] = new_tg

    add_unit_pair!(basis, pairset, s_ind, lm.deg)

    @assert s_ind == basis.input_load == ind_ord.max_ind

    return s_ind
end
