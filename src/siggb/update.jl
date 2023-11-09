# updating the pairset and basis

# add new reduced rows to basis/syzygies
# TODO: make a new symbol ht everytime?
function update_basis!(basis::Basis,
                       matrix::MacaulayMatrix,
                       pairset::Pairset{N},
                       symbol_ht::MonomialHashtable,
                       basis_ht::MonomialHashtable) where N

    new_basis_c = 0
    new_syz_c = 0

    # TODO: this can be done more optimally
    toadd = matrix.toadd[1:matrix.toadd_length]
    # sort!(toadd, by = i -> matrix.sigs[i], lt = lt_pot)

    @inbounds for i in toadd
        # determine if row is zero
        row = matrix.rows[i]
        new_sig = matrix.sigs[i]
        new_idx = index(new_sig)
        new_sig_mon = monomial(new_sig)
        new_sig_mask = (new_idx, divmask(new_sig_mon, basis_ht.divmap,
                                         basis_ht.ndivbits))
        if isempty(row)
            new_syz_c += 1

            # make sure we have enough space
            if basis.syz_load == basis.syz_size
                basis.syz_size *= 2
                resize!(basis.syz_sigs, basis.syz_size)
                resize!(basis.syz_masks, basis.syz_size)
            end
            
            # add new syz sig
            l = basis.syz_load + 1
            basis.syz_sigs[l] = monomial(new_sig)
            basis.syz_masks[l] = new_sig_mask
            basis.syz_load += 1

            # kill pairs with known syz signature
            @inbounds for j in 1:pairset.load
                p = pairset.elems[j]
                cond = index(p.top_sig) == new_idx
                if cond && divch(new_sig_mon, monomial(p.top_sig),
                                 new_sig_mask[2], p.top_sig_mask)
                    pairset.elems[j].top_index = 0
                end
                cond = index(p.bot_sig) == new_idx
                if cond && divch(new_sig_mon, monomial(p.bot_sig),
                                 new_sig_mask[2], p.bot_sig_mask)
                    pairset.elems[j].top_index = 0
                end
            end

            # remove pairs that became rewriteable in previous loop
            remove_red_pairs!(pairset)
        else
            new_basis_c += 1

            # make sure we have enough space
            if basis.basis_load == basis.basis_size
                basis.basis_size *= 2
                resize!(basis.sigs, basis.basis_size)
                resize!(basis.sigmasks, basis.basis_size)
                resize!(basis.sigratios, basis.basis_size)
                resize!(basis.rewrite_nodes, basis.basis_size)
                resize!(basis.lm_masks, basis.basis_size)
                resize!(basis.monomials, basis.basis_size)
                resize!(basis.coefficients, basis.basis_size)
                resize!(basis.is_red, basis.basis_size)
            end
            
            # add to basis hashtable
            insert_in_basis_hash_table_pivots!(row, basis_ht, symbol_ht)
            lm = basis_ht.exponents[first(row)]
            s = new_sig

            # add everything to basis
            l = basis.basis_load + 1
            basis.sigs[l] = new_sig
            basis.sigmasks[l] = new_sig_mask
            new_sig_ratio = divide(lm, new_sig_mon)
            basis.sigratios[l] = new_sig_ratio

            parent_ind = matrix.parent_inds[i]
            tree_data = basis.rewrite_nodes[parent_ind+1]
            insind = 3 
            @inbounds for j in insind:insind+tree_data[1]
                child_ind = tree_data[j]
                rat = basis.sigratios[child_ind-1]
                if lt_drl(new_sig_ratio, rat)
                    break
                end
                insind += 1
            end
            insert!(tree_data, insind, l+1)
            tree_data[1] += 1
            basis.rewrite_nodes[l+1] = [-1, parent_ind+1]
            
            basis.lm_masks[l] = divmask(lm, basis_ht.divmap,
                                        basis_ht.ndivbits)
            basis.monomials[l] = row
            basis.coefficients[l] = matrix.coeffs[i]
            basis.basis_load = l

            # build new pairs
            update_pairset!(pairset, basis, basis_ht, l)
        end
    end
    if new_basis_c != 0 || new_syz_c != 0
        @info "$(new_basis_c) new, $(new_syz_c) zero"
    end
end


# construct all pairs with basis element at new_basis_idx
# and perform corresponding rewrite checks
# assumes DPOT
function update_pairset!(pairset::Pairset{N},
                         basis::Basis,
                         basis_ht::MonomialHashtable,
                         new_basis_idx::Int) where N


    new_sig_mon = monomial(basis.sigs[new_basis_idx])
    new_sig_idx = index(basis.sigs[new_basis_idx])

    # check existing pairs for rewriteability against element
    # at new_basis_idx
    @inbounds bmask = basis.sigmasks[new_basis_idx]
    @inbounds parent_ind = basis.rewrite_nodes[new_basis_idx+1][2]
    @inbounds for i in 1:pairset.load
        p = pairset.elems[i]
        if p.top_index == parent_ind-1
            if divch(new_sig_mon, monomial(p.top_sig),
                     mask(bmask), p.top_sig_mask)
                pairset.elems[i].top_index = 0
                continue
            end
        end
        if p.bot_index == parent_ind-1
            if divch(new_sig_mon, monomial(p.bot_sig),
                     mask(bmask), p.bot_sig_mask)
                pairset.elems[i].top_index = 0
            end
        end
    end

    # kill pairs that became rewriteable in the previous round
    remove_red_pairs!(pairset)

    # resize pairset if needed
    num_new_pairs = new_basis_idx - 1
    if pairset.load + num_new_pairs >= pairset.size
          resize!(pairset.elems, max(2 * pairset.size,
                                     pairset.load - num_new_pairs))
          pairset.size *= 2
    end

    new_sig_ratio = basis.sigratios[new_basis_idx]
    new_lm = leading_monomial(basis, basis_ht, new_basis_idx)
    # pair construction loop
    @inbounds for i in basis.basis_offset:(new_basis_idx - 1)
        basis_lm = leading_monomial(basis, basis_ht, i)
        basis_sig_idx = index(basis.sigs[i])

        mult_new_elem = lcm_div(new_lm, basis_lm)
        new_pair_sig_mon = mul(mult_new_elem, new_sig_mon)
        
        mult_basis_elem = lcm_div(basis_lm, new_lm)
        basis_pair_sig_mon = mul(mult_basis_elem, monomial(basis.sigs[i]))

        new_pair_sig = (new_sig_idx, new_pair_sig_mon)
        basis_pair_sig = (basis_sig_idx, basis_pair_sig_mon)
        
        # check if S-pair is singular
        new_pair_sig == basis_pair_sig && continue

        new_pair_sig_mask = divmask(new_pair_sig_mon,
                                    basis_ht.divmap,
                                    basis_ht.ndivbits)
        basis_pair_sig_mask = divmask(basis_pair_sig_mon,
                                      basis_ht.divmap,
                                      basis_ht.ndivbits)

        # check both pair sigs against non-trivial syzygies
        rewriteable_syz(basis, new_pair_sig,
                        new_pair_sig_mask) && continue
        rewriteable_syz(basis, basis_pair_sig,
                        basis_pair_sig_mask) && continue

        # check both pair signatures against koszul syzygies
        rewriteable_koszul(basis, basis_ht, new_pair_sig,
                           new_pair_sig_mask) && continue
        rewriteable_koszul(basis, basis_ht, basis_pair_sig,
                           basis_pair_sig_mask) && continue

        top_sig, top_sig_mask, top_index,
        bot_sig, bot_sig_mask, bot_index = begin
	    if lt_pot(basis_pair_sig, new_pair_sig)
                new_pair_sig, new_pair_sig_mask, new_basis_idx,
                basis_pair_sig, basis_pair_sig_mask, i
            else
                basis_pair_sig, basis_pair_sig_mask, i,
                new_pair_sig, new_pair_sig_mask, new_basis_idx
            end 
        end
        
        pair_deg = new_pair_sig_mon.deg + basis.degs[new_sig_idx]
        new_pair =  SPair(top_sig, bot_sig,
                          top_sig_mask, bot_sig_mask,
                          top_index, bot_index,
                          pair_deg)
        pairset.elems[pairset.load + 1] = new_pair
        pairset.load += 1
    end
end

@inline function rewriteable_syz(basis::Basis,
                                 sig::Sig,
                                 sigmask::DivMask)

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
                                   sigmask::DivMask)

    k = find_canonical_rewriter(basis, sig, sigmask)
    return k != idx
end

# @inline function rewriteable_basis(basis::Basis,
#                                    idx::Int,
#                                    sig::Sig,
#                                    sigmask::DivMask)

#     ind = index(sig)
    
#     @inbounds for i in basis.basis_load:-1:basis.basis_offset
#         i == idx && continue
#         i_sig_idx = index(basis.sigmasks[i])
#         i_sig_idx != ind && continue
#         i_sig_mask = mask(basis.sigmasks[i])
#         if divch(monomial(basis.sigs[i]), monomial(sig),
#                  i_sig_mask, sigmask)
#             is_rewr = comp_sigratio(basis, i, idx)
#             is_rewr && return true
#         end
#     end
#     return false
# end

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

@inline function rewriteable_koszul(basis::Basis,
                                    basis_ht::MonomialHashtable,
                                    sig::Sig,
                                    sigmask::DivMask)

    @inbounds for i in basis.basis_offset:basis.basis_load
        if index(basis.sigs[i]) < index(sig)
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
                     sigmask::DivMask)

    rewriteable_syz(basis, sig, sigmask) && return true
    rewriteable_basis(basis, idx, sig, sigmask) && return true
    rewriteable_koszul(basis, basis_ht, sig, sigmask) && return true
    return false
end

# helper functions for readability
function leading_monomial(basis::Basis,
                          basis_ht::MonomialHashtable,
                          i)

    return basis_ht.exponents[first(basis.monomials[i])]
end

@inline function comp_sigratio(basis::Basis, ind1::Int, ind2::Int)

    rat1 = basis.sigratios[ind1]
    rat2 = basis.sigratios[ind2]
    if rat1 == rat2
        return lt_drl(monomial(basis.sigs[ind1]), monomial(basis.sigs[ind2]))
    end
    return lt_drl(rat1, rat2)
end

# remove pairs that are rewriteable
function remove_red_pairs!(pairset::Pairset)
    iszero(pairset.load) && return
    j = 0 
    @inbounds for i in 1:pairset.load
        iszero(pairset.elems[i].top_index) && continue
        j += 1
        pairset.elems[j] = pairset.elems[i]
    end
    pairset.load = j 
end
