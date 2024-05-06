function construct_module_wrap(sig::Sig{N},
                               basis::Basis{N},
                               basis_ht::MonomialHashtable{N},
                               mat_index::Int,
                               tr::Tracer,
                               vchar::Val{Char},
                               ind_ord::IndOrder,
                               idx::SigIndex) where {N, Char}

    cofac_coeffs, cofac_mons_hsh = construct_module(sig, basis,
                                                    basis_ht,
                                                    mat_index,
                                                    tr, vchar,
                                                    ind_ord, idx)

    if !isempty(cofac_coeffs)
        sort_poly!((cofac_coeffs, cofac_mons_hsh),
                   by = midx -> basis_ht.exponents[midx],
                   lt = lt_drl, rev = true)
        # normalize cofac coefficients
        normalize_cfs!(cofac_coeffs, vchar)
    end

    return cofac_coeffs, cofac_mons_hsh
end

    
function construct_module(basis::Basis{N},
                          basis_ht::MonomialHashtable{N},
                          basis_index::Int,
                          tr::Tracer,
                          vchar::Val{Char},
                          ind_order::IndOrder,
                          idx::SigIndex) where {N, Char}

    @inbounds sig = basis.sigs[basis_index]

    if basis_index >= basis.basis_offset
        ik = basis.mod_rep_known[basis_index]
        if ik[idx]
            return basis.mod_reps[basis_index][idx]
        end

        @inbounds mat_ind = tr.basis_ind_to_mat[basis_index]
        res = construct_module_core(sig, basis, basis_ht,
                                    mat_ind, tr,
                                    vchar, 
                                    ind_order, idx)

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
                          idx::SigIndex) where {N, Char}
    
    if ind_ord.ord[index(sig)] < ind_ord.ord[idx]
        return Coeff[], MonIdx[]
    end

    tr_mat = tr.mats[mat_index]

    row_ind, rewr_basis_ind = tr_mat.rows[sig]

    basis_ind = get(tr_mat.is_basis_row, row_ind, 0)
    if !iszero(basis_ind)
        @assert basis.sigs[basis_ind] == sig
        return construct_module(basis, basis_ht, basis_ind,
                                tr, vchar, ind_ord, idx)
    end

    return construct_module_core(sig, basis, basis_ht, mat_index,
                                 tr, vchar, ind_ord, idx)
end

function construct_module_core(sig::Sig{N},
                               basis::Basis{N},
                               basis_ht::MonomialHashtable{N},
                               mat_index::Int,
                               tr::Tracer,
                               vchar::Val{Char},
                               ind_ord::IndOrder,
                               idx::SigIndex) where {N, Char}

    if ind_ord.ord[index(sig)] < ind_ord.ord[idx]
        return Coeff[], MonIdx[]
    end

    tr_mat = tr.mats[mat_index]

    row_ind, rewr_basis_ind = tr_mat.rows[sig]

    # construct module representation of canonical rewriter
    rewr_mod_cfs, rewr_mod_mns = construct_module(basis, basis_ht,
                                                  rewr_basis_ind,
                                                  tr, vchar,
                                                  ind_ord, idx)

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
        sort_poly!((res_mod_cfs, res_mod_mns))
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
                                     idx)
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
