function construct_module(basis::Basis{N},
                          basis_index::Int,
                          tr::Tracer,
                          vchar::Val{Char},
                          mod_cache::Dict{Sig, Vector{Polynomial{N}}},
                          mod_dim::SigIndex,
                          ind_order::IndOrder,
                          just_index::SigIndex=SigIndex(0)) where {N, Char}

    @inbounds sig = basis.sigs[basis_index]

    if basis_index >= basis.basis_offset
        @inbounds mat_ind = tr.basis_ind_to_mat[basis_index]
        return construct_module(sig, basis, mat_ind, tr,
                                vchar, mod_cache, 
                                mod_dim, ind_order,
                                just_index)
    else
        # if it was an input element we just take the signature
        res = [(Coeff[], Monomial{N}[]) for _ in 1:mod_dim]
        if iszero(just_index) || index(sig) == just_index
            res[index(sig)] = ([one(Coeff)], [monomial(sig)])
        end
        return res
    end
        
end

# construct a module representation out of a given sig
function construct_module(sig::Sig{N},
                          basis::Basis{N},
                          mat_index::Int,
                          tr::Tracer,
                          vchar::Val{Char},
                          mod_cache::Dict{Sig, Vector{Polynomial{N}}},
                          mod_dim::SigIndex,
                          ind_ord::IndOrder,
                          just_index::SigIndex=SigIndex(0)) where {N, Char}

    if haskey(mod_cache, sig)
        return mod_cache[sig]
    end

    if ind_ord.ord[index(sig)] < ind_ord.ord[just_index]
        return [(Coeff[], Monomial{N}[]) for _ in 1:mod_dim]
    end

    tr_mat = tr.mats[mat_index]

    row_ind, rewr_basis_ind = tr_mat.rows[sig]

    # construct module representation of canonical rewriter
    rewr_mod = construct_module(basis, rewr_basis_ind, tr, vchar,
                                mod_cache,
                                mod_dim, ind_ord, just_index)

    # multiply by monomial
    rewr_sig = basis.sigs[rewr_basis_ind]
    mult = divide(monomial(sig), monomial(rewr_sig))
    isone = all(iszero, mult.exps)
    res = Vector{Polynomial{N}}(undef, mod_dim)
    @inbounds for i in 1:mod_dim
        res[i] = (copy(rewr_mod[i][1]), isone ? copy(rewr_mod[i][2]) : mul_by_mon(rewr_mod[i][2], mult))
    end

    @inbounds row_ops = tr_mat.col_inds_and_coeffs[row_ind]
    @inbounds for (j, coeff)  in row_ops
        j_sig = tr_mat.row_ind_to_sig[j]
        j_sig_mod = construct_module(j_sig, basis, mat_index,
                                     tr, vchar, mod_cache,
                                     mod_dim, ind_ord, just_index)
        for i in 1:mod_dim
            !iszero(just_index) && i != just_index
            res_i_coeffs = res[i][1]
            res_i_mons = res[i][2]
            j_mod_coeffs = j_sig_mod[i][1]
            mul_j_mod_coeffs = mul_by_coeff(j_mod_coeffs, addinv(coeff, vchar), vchar)
            j_mod_mons = j_sig_mod[i][2]
            res[i] = add_pols(res_i_mons, j_mod_mons,
                              res_i_coeffs, mul_j_mod_coeffs,
                              vchar)
        end
    end

    diag_coeff = tr_mat.diagonal[row_ind]
    @inbounds for i in 1:mod_dim
        res_i_coeffs = res[i][1]
        mul_by_coeff!(res_i_coeffs, diag_coeff, vchar)
    end

    if haskey(mod_cache, sig)
        @assert mod_cache[sig] == res
    end
    mod_cache[sig] = res
    return res
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

function add_pols(mons1::Vector{M},
                  mons2::Vector{M},
                  coeffs1::Vector{Coeff},
                  coeffs2::Vector{Coeff},
                  vch::Val{Char}) where {M <: Monomial, Char}

    l1 = length(mons1)
    l2 = length(mons2)
    mons_res = Vector{M}(undef, l1 + l2)
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
        elseif lt_drl(m2, m1)
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

    return coeffs_res, mons_res
end
