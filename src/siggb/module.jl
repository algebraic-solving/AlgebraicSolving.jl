# construct a module representation out of a given sig
function construct_module(sig::Sig{N},
                          basis::Basis{N},
                          tr::Tracer,
                          vchar::Val{Char},
                          R::MPolyRing) where {N, Char}

    degs = basis.degs
    mod_dim = length(degs)

    @inbounds deg = monomial(sig).deg + degs[index(sig)]
    tr_mat = tr[deg]

    row_ind, rewr_basis_ind = tr_mat.rows[sig]
    if rewr_basis_ind >= basis.basis_offset
        @inbounds rewr_sig = basis.sigs[rewr_basis_ind]

        # construct module representation of canonical rewriter
        res = construct_module(rewr_sig, basis, tr, vchar, R)

        # multiply by monomial
        mult = divide(monomial(sig), monomial(rewr_sig))
        @inbounds for i in 1:mod_dim
            res_i_mons = res[i][2]
            mul_by_mon!(res_i_mons, mult)
        end
    else
        # if it was an input element we just take the signature
        res = [(Coeff[], Monomial{N}[]) for _ in 1:mod_dim]
        res[index(sig)] = ([one(Coeff)], [monomial(sig)])
    end

    @inbounds row_ops = tr_mat.col_inds_and_coeffs[row_ind]
    @inbounds for (j, coeff)  in row_ops
        j_sig = tr_mat.row_ind_to_sig[j]
        j_sig_mod = construct_module(j_sig, basis, tr, vchar, R)
        for i in 1:mod_dim
            res_i_coeffs = res[i][1]
            res_i_mons = res[i][2]
            j_mod_coeffs = j_sig_mod[i][1]
            mul_by_coeff!(j_mod_coeffs, addinv(coeff, vchar), vchar)
            j_mod_mons = j_sig_mod[i][2]
            res[i] = add_pols(res_i_mons, j_mod_mons,
                              res_i_coeffs, j_mod_coeffs,
                              vchar, R)
        end
    end

    diag_coeff = tr_mat.diagonal[row_ind]
    @inbounds for i in 1:mod_dim
        res_i_coeffs = res[i][1]
        mul_by_coeff!(res_i_coeffs, diag_coeff, vchar)
    end

    return res
end

# functions for polynomials
function mul_by_mon!(mons::Vector{M},
                     mon::M) where {M <: Monomial}

    @inbounds for i in 1:length(mons)
        mons[i] = mul(mon, mons[i])
    end
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
                  vch::Val{Char},
                  R::MPolyRing) where {M <: Monomial, Char}

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
