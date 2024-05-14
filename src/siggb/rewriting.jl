@inline function rewriteable_syz(basis::Basis,
                                 sig::Sig,
                                 sigmask::DivMask,
                                 tags::Tags)

    ind = index(sig)

    @inbounds for i in 1:basis.syz_load
        if index(basis.syz_masks[i]) == ind
            is_rewr = divch(basis.syz_sigs[i], monomial(sig),
                            basis.syz_masks[i][2], sigmask) 
            is_rewr && return true
        end
    end
    return false
end

@inline function rewriteable_basis(basis::Basis,
                                   idx::Int,
                                   sig::Sig,
                                   sigmask::DivMask,
                                   tags::Tags)

    ind = index(sig)

    k = find_canonical_rewriter(basis, sig, sigmask)
    return k != idx
end

function find_canonical_rewriter(basis::Basis,
                                 sig::Sig,
                                 sigmask::DivMask)

    node_ind = 1
    while true
        @inbounds node = basis.rewrite_nodes[node_ind]
        @inbounds node[1] == -1 && break
        found_div = false
        @inbounds for i in 3:3+node[1]
            ch = node[i]
            basis_sig = basis.sigs[ch - 1]
            basis_sigmask = basis.sigmasks[ch - 1]
            index(basis_sig) != index(sig) && continue
            if divch(monomial(basis_sig), monomial(sig),
                     mask(basis_sigmask), sigmask)
                node_ind = ch
                found_div = true
                break
            end
        end
        !found_div && break
    end
    return node_ind - 1
end

# function find_canonical_rewriter(basis::Basis,
#                                  sig::Sig,
#                                  sigmask::DivMask)

#     cand = zero(Int)
#     for i in 1:basis.basis_load
#         basis.input_load < i < basis.basis_offset && continue
#         basis.is_red[i] && continue
#         bs = basis.sigs[i]
#         if index(bs) == index(sig) && divch(monomial(bs), monomial(sig),
#                                             mask(basis.sigmasks[i]), sigmask)
#             if iszero(cand) || comp_sigratio(basis, i, cand)
#                 cand = i
#             end
#         end
#     end
#     return cand
# end

@inline function rewriteable_koszul(basis::Basis,
                                    basis_ht::MonomialHashtable,
                                    sig::Sig,
                                    sigmask::DivMask,
                                    ind_order::IndOrder,
                                    tags::Tags)

    s_ind = index(sig)

    @inbounds for i in basis.basis_offset:basis.basis_load
        b_ind = index(basis.sigs[i])
        if (ind_order.ord[b_ind] < ind_order.ord[s_ind]
            && !are_incompat(b_ind, s_ind, ind_order))
            if divch(basis.lm_masks[i], sigmask)
                if divch(leading_monomial(basis, basis_ht, i), monomial(sig))
                    return true
                end
            end
        end
    end
    return false
end

function rewriteable(basis::Basis,
                     basis_ht::MonomialHashtable,
                     idx::Int,
                     sig::Sig,
                     sigmask::DivMask,
                     ind_order::IndOrder,
                     tags::Tags)

    s_ind = index(sig)

    rewriteable_syz(basis, sig, sigmask, tags) && return true
    rewriteable_basis(basis, idx, sig, sigmask, tags) && return true
    rewriteable_koszul(basis, basis_ht, sig, sigmask, ind_order, tags) && return true
    return false
end

@inline function comp_sigratio(basis::Basis, ind1::Int, ind2::Int)

    rat1 = basis.sigratios[ind1]
    rat2 = basis.sigratios[ind2]
    if rat1 == rat2
        return true
    end
    return lt_drl(rat1, rat2)
end
