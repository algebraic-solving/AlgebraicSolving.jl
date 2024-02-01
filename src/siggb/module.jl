function construct_module(basis::Basis{N},
                          basis_ht::MonomialHashtable{N},
                          basis_index::Int,
                          tr::Tracer,
                          vchar::Val{Char},
                          ind_order::IndOrder,
                          idx::SigIndex,
                          gb::Vector{T};
                          maintain_nf::Bool=false) where {N, Char,
                                                          T <: MPolyRingElem}

    @inbounds sig = basis.sigs[basis_index]

    if basis_index >= basis.basis_offset
        if isassigned(basis.mod_rep_known, basis_index) 
            ik = basis.mod_rep_known[basis_index]
            if ik[idx]
                return basis.mod_reps[basis_index][idx]
            end
        else
            basis.mod_rep_known[basis_index] = falses(ind_order.max_ind)
            basis.mod_reps[basis_index] =
                Vector{Polynomial}(undef, ind_order.max_ind)
        end

        @inbounds mat_ind = tr.basis_ind_to_mat[basis_index]
        res = construct_module_core(sig, basis, basis_ht,
                                    mat_ind, tr,
                                    vchar, 
                                    ind_order, idx,
                                    gb,
                                    maintain_nf=maintain_nf)

        if maintain_nf
            res_pol = convert_to_pol(parent(first(gb)),
                                     [basis_ht.exponents[midx] for midx in res[2]],
                                     res[1])
            f = open("./nf_log.txt", "a+")
            println(f, res_pol)
            close(f)
            res_pol_nf = my_normal_form([res_pol], gb)[1]
            res = convert_to_ht(res_pol_nf, basis_ht, vchar, normalise=false)
        end
        basis.mod_rep_known[basis_index][idx] = true
        basis.mod_reps[basis_index][idx] = res
        return res
    else
        # if it was an input element we just take the signature
        if index(sig) == idx
            check_enlarge_hashtable!(basis_ht, 1)
            mon_hsh_idx = insert_in_hash_table!(basis_ht, monomial(sig))
            return ([one(Coeff)], [mon_hsh_idx])
        else
            return (Coeff[], MonIdx[])
        end
    end
end

function construct_module(sig::Sig{N},
                          basis::Basis{N},
                          basis_ht::MonomialHashtable{N},
                          mat_index::Int,
                          tr::Tracer,
                          vchar::Val{Char},
                          ind_ord::IndOrder,
                          idx::SigIndex,
                          gb::Vector{T};
                          maintain_nf::Bool=false) where {N, Char, T <: MPolyRingElem}
    
    if ind_ord.ord[index(sig)] < idx
        return Coeff[], Monomial{N}[]
    end

    tr_mat = tr.mats[mat_index]

    row_ind, rewr_basis_ind = tr_mat.rows[sig]

    basis_ind = get(tr_mat.is_basis_row, row_ind, 0)
    if !iszero(basis_ind)
        @assert basis.sigs[basis_ind] == sig
        return construct_module(basis, basis_ht, basis_ind,
                                tr, vchar, ind_ord, idx, gb,
                                maintain_nf = maintain_nf)
    end

    return construct_module_core(sig, basis, basis_ht, mat_index,
                                 tr, vchar, ind_ord, idx, gb,
                                 maintain_nf = maintain_nf)
end

function construct_module_core(sig::Sig{N},
                               basis::Basis{N},
                               basis_ht::MonomialHashtable{N},
                               mat_index::Int,
                               tr::Tracer,
                               vchar::Val{Char},
                               ind_ord::IndOrder,
                               idx::SigIndex,
                               gb::Vector{T};
                               maintain_nf::Bool=false) where {N, Char,
                                                               T <: MPolyRingElem}

    if ind_ord.ord[index(sig)] < idx
        return Coeff[], Monomial{N}[]
    end

    tr_mat = tr.mats[mat_index]

    row_ind, rewr_basis_ind = tr_mat.rows[sig]

    # construct module representation of canonical rewriter
    rewr_mod_cfs, rewr_mod_mns = construct_module(basis, basis_ht,
                                                  rewr_basis_ind,
                                                  tr, vchar,
                                                  ind_ord, idx, gb,
                                                  maintain_nf=maintain_nf)

    # multiply by monomial
    rewr_sig = basis.sigs[rewr_basis_ind]
    mult = divide(monomial(sig), monomial(rewr_sig))
    plength = length(rewr_mod_cfs)
    res_mod_cfs, res_mod_mns = copy(rewr_mod_cfs), Vector{MonIdx}(undef, plength)
    if all(iszero, mult.exps)
        @inbounds for i in 1:plength
            res_mod_mns[i] = rewr_mod_mns[i]
        end
    else
        hsh = Base.hash(mult)
        check_enlarge_hashtable!(basis_ht, plength)
        insert_multiplied_poly_in_hash_table!(res_mod_mns, hsh, mult,
                                              rewr_mod_mns,
                                              basis_ht, basis_ht)
        s = sortperm(res_mod_mns)
        res_mod_cfs = res_mod_cfs[s]
        res_mod_mns = res_mod_mns[s]
    end

    # construct module rep of all reducers
    @inbounds row_ops = tr_mat.col_inds_and_coeffs[row_ind]
    @inbounds for (j, coeff)  in row_ops
        j_sig = tr_mat.row_ind_to_sig[j]
        if ind_ord.ord[index(j_sig)] < ind_ord.ord[idx]
            continue
        end
        j_sig_mod = construct_module(j_sig, basis, basis_ht,
                                     mat_index,
                                     tr, vchar,
                                     ind_ord,
                                     idx, gb,
                                     maintain_nf=maintain_nf)
        j_mod_coeffs = j_sig_mod[1]
        mul_j_mod_coeffs = mul_by_coeff(j_mod_coeffs, addinv(coeff, vchar),
                                        vchar)
        j_mod_mons = j_sig_mod[2]
        res_mod_cfs, res_mod_mns = add_pols(res_mod_cfs, res_mod_mns,
                                            mul_j_mod_coeffs, j_mod_mons,
                                            vchar)
    end

    diag_coeff = tr_mat.diagonal[row_ind]
    mul_by_coeff!(res_mod_cfs, diag_coeff, vchar)

    return res_mod_cfs, res_mod_mns
end

# functions for polynomials
function mul_by_mon(mons::Vector{M},
                    mon::M) where {M <: Monomial}

    mons_res = Vector{M}(undef, length(mons))
    @inbounds for i in 1:length(mons)
        mons_res[i] = mul(mon, mons[i])
    end
    return mons_res
end

function mul_by_coeff(coeffs::Vector{Coeff},
                      c::Coeff,
                      vchar::Val{Char}) where Char 

    coeffs_res = Vector{Coeff}(undef, length(coeffs))
    @inbounds for i in 1:length(coeffs)
        coeffs_res[i] = mul(c, coeffs[i], vchar)
    end
    return coeffs_res
end

function mul_by_coeff!(coeffs::Vector{Coeff},
                       c::Coeff,
                       vchar::Val{Char}) where Char 

    @inbounds for i in 1:length(coeffs)
        coeffs[i] = mul(c, coeffs[i], vchar)
    end
end

# assumes mons are sorted ascendingly
function add_pols(coeffs1::Vector{Coeff},
                  mons1::Vector{MonIdx},
                  coeffs2::Vector{Coeff},
                  mons2::Vector{MonIdx},
                  vch::Val{Char}) where {Char}

    l1 = length(mons1)
    l2 = length(mons2)
    mons_res = Vector{MonIdx}(undef, l1 + l2)
    coeffs_res = Vector{Coeff}(undef, l1 + l2)
    
    ind1 = 1
    ind2 = 1
    new_l = 0 
    @inbounds while ind1 <= l1 && ind2 <= l2
        new_l += 1
        m1 = mons1[ind1]
        m2 = mons2[ind2]
        if m1 == m2
            mons_res[new_l] = m1
            coeffs_res[new_l] = add(coeffs1[ind1], coeffs2[ind2], vch)
            ind1 += 1
            ind2 += 1
        elseif m1 < m2
            mons_res[new_l] = m1
            coeffs_res[new_l] = coeffs1[ind1]
            ind1 += 1
        else
            mons_res[new_l] = m2
            coeffs_res[new_l] = coeffs2[ind2]
            ind2 += 1
        end
    end

    while ind1 <= l1
        new_l += 1
        mons_res[new_l] = mons1[ind1]
        coeffs_res[new_l] = coeffs1[ind1]
        ind1 += 1
    end

    while ind2 <= l2 
        new_l += 1
        mons_res[new_l] = mons2[ind2]
        coeffs_res[new_l] = coeffs2[ind2]
        ind2 += 1
    end

    resize!(mons_res, new_l)
    resize!(coeffs_res, new_l)

    zero_cfs_inds = findall(c -> iszero(c), coeffs_res)
    deleteat!(mons_res, zero_cfs_inds)
    deleteat!(coeffs_res, zero_cfs_inds)

    return coeffs_res, mons_res
end


function add_pols_2(cfs1::Vector{Coeff},
                    exps1::Vector{Monomial{N}},
                    cfs2::Vector{Coeff},
                    exps2::Vector{Monomial{N}},
                    char::Val{Char}) where {N, Char}

    R, vrs = polynomial_ring(GF(Int(Char)), ["x$i" for i in 1:N],
                             ordering = :degrevlex)
    p1 = convert_to_pol(R, exps1, cfs1)
    p2 = convert_to_pol(R, exps2, cfs2)
    p = p1+p2

    lp = length(p)
    exps = exponent_vectors(p)
    cfs = coefficients(p)
    
    res_exps = Vector{Monomial{N}}(undef, lp)
    res_cfs = Vector{Coeff}(undef, lp)
    @inbounds for (i, (cf, evec)) in enumerate(zip(cfs, exps)) 
        m = monomial(SVector{N}((Exp).(evec)))
        cff = Int(lift(ZZ, cf))
        res_exps[i] = m
        res_cfs[i] = cff
    end

    return res_cfs, res_exps
end
