# updating the pairset and basis

# for readability
function leading_monomial(basis::POTBasis,
                          basis_ht::MonomialHashtable,
                          i)

    return basis_ht.exponents[first(basis.monomials[i])]
end

# construct all pairs with basis element at new_basis_idx
# and perform corresponding rewrite checks
# assumes DPOT
function update_pairset!(pairset::Pairset{SPair{N}},
                         basis::Basis,
                         basis_ht::MonomialHashtable,
                         new_basis_idx::Int) where N


    new_sig_mon = monomial(basis.sigs[new_basis_idx])
    new_sig_idx = index(basis.sigs[new_basis_idx])

    # check existing pairs for rewriteability against element
    # at new_basis_idx
    bmask = basis.sigmasks[new_basis_idx]
    @inbounds for i in 1:pairset.load
        p = pairset.pairs[i]
        (iszero(p.top_index) || index(p.top_sig) != new_sig_idx) && continue
        if div(new_sig_mon, p.top_sig, bmask, p.top_sig_mask)
            pairset.pairs[i].top_index = 0
            continue
        end
        new_sig_idx != index(p.bot_sig) && continue
        if div(new_sig_mon, p.bot_sig, bmask, p.bot_sig_mask)
            pairset.pairs[i].top_index = 0
            continue
        end
    end

    # remove pairs that are rewriteable
    j = 1
    @inbounds for i in 1:pairset.load
        iszero(pairset.pairs[i].top_index) && continue
        pairset[j] = pairset[i]
        j += 1
    end
    pairset.load -= (pairset.load - j)

    # resize pairset if needed
    pair_size = length(pairset.pairs)
    num_new_pairs = new_basis_idx - 1
    if pairset.load + num_new_pairs >= pair_size
        resize!(pairset.pairs, max(2 * pair_size, pair.load - num_new_pairs))
    end

    new_lm = leading_monomial(basis, basis_ht, new_basis_idx)
    # pair construction loop
    @inbounds for i in 1:(new_basis_idx - 1)
        basis_lm = basis_ht.exponents[basis.lms[i]]
        ind = index(basis.sigs[i])

        mult_new_elem = lcm_div(new_lm, basis_lm)
        new_pair_sig_mon = mul(mult_new_elem, new_sig_mon)

        mult_basis_elem = lcm_div(basis_lm, new_lm)
        basis_pair_sig_mon = mul(mult_basis_elem, monomial(basis.sigs[i]))

        # check if S-pair is singular
        new_sig_idx == ind && new_pair_sig_mon == basis_pair_sig_mon && continue

        is_rewr = false
        new_pair_sig_mask = divmask(new_pair_sig_mon,
                                    basis_ht.divmap,
                                    basis_ht.ndivbits)
        basis_pair_sig_mask = divmask(basis_pair_sig_mon,
                                      basis_ht.divmap,
                                      basis_ht.ndivbits)

        # check both pair sigs against non-trivial syzygies
        @inbounds for j in 1:basis.syz_load
            if index(basis.syz_masks[j]) == new_sig_idx
                is_rewr = div(basis.syz_sigs[j], new_pair_sig_mon,
                              basis.syz_masks[j][2], new_pair_sig_mask) 
                is_rewr && break
            end
            if index(basis.syz_masks[j]) == ind
                is_rewr = div(basis.syz_sigs[j], basis_pair_sig_mon,
                              basis.syz_masks[j][2], basis_pair_sig_mask) 
                is_rewr && break
            end
        end
        is_rewr && continue

        # check both pair signatures against koszul syzygies
        # TODO: should we store the indices with the lm masks
        @inbounds for j in 1:(new_basis_idx-1)
            if index(basis.sigs[j]) < new_basis_idx
                is_rewr = div(leading_monomial(basis, basis_ht, j),
                              new_pair_sig_mon,
                              basis.lm_masks[j], new_pair_sig_mask)
                is_rewr && break
            end
            if index(basis.sigs[j]) < ind
                is_rewr = div(leading_monomial(basis, basis_ht, j),
                              basis_pair_sig_mon,
                              basis.lm_masks[j], basis_pair_sig_mask)
                is_rewr && break
            end
        end
        is_rewr && continue
        
        # check multiplied signature of basis element against basis sigs
        @inbounds for j in 1:new_basis_idx
            index(basis.sigmasks[j]) != new_sig_idx && continue
            is_rewr = div(monomial(basis.sigs[j]), basis_pair_sig_mon,
                          basis.sigmasks[j][2], basis_pair_sig_mask)
            is_rewr && break
        end
        is_rewr && continue
        
        # TODO: feels like this could be simplified
        new_pair = if index(basis.sigs[i]) == new_sig_idx
            if lt_drl(new_sig_mon, basis_pair_sig_mon)
                SPair(new_sig_mon, Sig(ind, basis_pair_sig_mon),
                      new_pair_sig_mask, basis_pair_sig_mask, new_basis_idx, i)
            else
                SPair(basis_pair_sig_mon, new_sig_mon,
                      basis_pair_mon_mask, new_pair_sig_mask, i, new_basis_idx)
            end
        elseif index(basis.sigs[i]) > new_sig_idx
            SPair(basis_pair_sig_mon, new_sig_mon,
                     basis_pair_mon_mask, new_pair_sig_mask, i, new_basis_idx)
        else
            SPair(new_sig_mon, Sig(ind, basis_pair_sig_mon),
                  new_pair_sig_mask, basis_pair_sig_mask, new_basis_idx, i)
        end
            
        pairset.pairs[pairset.load + 1] = new_pair
        pairset.load += 1
    end
end
