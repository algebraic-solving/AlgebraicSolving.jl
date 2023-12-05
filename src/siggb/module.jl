# construct a module representation out of a given sig
function construct_module(sig::Sig{N},
                          degs::Vector{Exp},
                          basis::Basis{N},
                          tr::Tracer) where N

    mod_dim = length(degs)

    res = [(Coeff[], Monomial{N}[]) for _ in 1:mod_dim]

    deg, row_ind = find_deg_and_row_ind(sig, tr, degs)

    # construct module representation of canonical rewriter
    # TODO: check that bas_ind[row_ind] is not zero
    @inbounds rewr_sig = basis.sigs[bas_ind[row_ind]]
    mod_rep_rewr = construct_module(rewr_sig, degs, basis, tr)

    # TODO: multiply by monomial and coefficient
    @inbounds for i in 1:mod_dim
        res[i] = mod_rep_rewr[i]
    end

    # check what row reductions we did
    tr_mat = tr.matrices[deg]
    tr_mat_ind = get(tr_mat.row_inds, row_ind, 0)

    if iszero(tr_mat_ind)
        return res
    end

    @inbounds col_inds_coeffs = tr_mat.col_inds_and_coeffs[tr_mat_ind]
    @inbounds for (j, coeff)  in col_inds_coeffs
        j_sig = tr.sigs[deg][j]
        j_sig_mod = construct_module(sig, degs, basis, tr)
        # TODO: add coeff*j_sig_mod to res
        # TODO: are we sure that everything is sorted
    end

    return res
end

function find_deg_and_row_ind(sig::Sig,
                              tr::Tracer,
                              degs::Vector{Exp})

    @inbounds deg = degree(monomial(sig)) + degs[index(sig)]
    tr_ind = deg

    row_ind = 0
    @inbounds for (i, s) in enumerate(tr.sigs[deg])
        if s == sig
            row_ind = i
            break
        end
    end

    iszero(row_ind) && error("sig not found in matrix")

    return deg, row_ind
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
                       vchar::Val{Char}) 

    @inbounds for i in 1:length(coeffs)
        coeffs[i] = mul(c, coeffs[i], varch)
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

    return mons_res, coeffs_res
end
