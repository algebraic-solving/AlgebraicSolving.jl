# updating the pairset and basis

# for readability
function leading_monomial(basis::POTBasis,
                          basis_ht::MonomialHashtable,
                          i)

    return basis_ht.exponents[first(basis.monomials[i])]
end

@inline function comp_sigratio(basis::Basis, ind1::Int, ind2::Int2)

    rat1 = basis.sigratios[ind1]
    rat2 = basis.sigratios[ind2]
    if rat1 == rat2
        return lt_drl(monomial(basis.sigs[ind1]), monomial(basis.sigs[ind2]))
    end
    return lt_drl(rat1, rat2)
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
            if comp_sigratio(basis, new_sig_idx, p.top_index)
                pairset.pairs[i].top_index = 0
                continue
            end
        end
        new_sig_idx != index(p.bot_sig) && continue
        if div(new_sig_mon, p.bot_sig, bmask, p.bot_sig_mask)
            if comp_sigratio(basis, new_sig_idx, p.bot_index)
                pairset.pairs[i].top_index = 0
                continue
            end
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

    new_sig_ratio = basis.sigratios[new_basis_idx]
    new_lm = leading_monomial(basis, basis_ht, new_basis_idx)
    # pair construction loop
    @inbounds for i in 1:(new_basis_idx - 1)
        basis_lm = basis_ht.exponents[basis.lms[i]]
        basis_sig_idx = index(basis.sigs[i])

        mult_new_elem = lcm_div(new_lm, basis_lm)
        new_pair_sig_mon = mul(mult_new_elem, new_sig_mon)

        mult_basis_elem = lcm_div(basis_lm, new_lm)
        basis_pair_sig_mon = mul(mult_basis_elem, monomial(basis.sigs[i]))

        # check if S-pair is singular
        new_sig_idx == basis_sig_idx && new_pair_sig_mon == basis_pair_sig_mon && continue

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
            if index(basis.syz_masks[j]) == basis_sig_idx
                is_rewr = div(basis.syz_sigs[j], basis_pair_sig_mon,
                              basis.syz_masks[j][2], basis_pair_sig_mask) 
                is_rewr && break
            end
        end
        is_rewr && continue

        # check both pair sigs against basis sigs
        @inbounds for j in 1:basis.basis_load-1
            j == new_basis_idx && continue
            j_sig_idx = index(basis.sigmasks[j])
            (j_sig_idx != new_sig_idx && j_sig_idx != basis_sig_idx) && continue
            j_sig_mask = basis.sigmasks[j][2]
            if (j_sig_idx == new_sig_idx &&
                div(monomial(basis.sigs[j]), new_pair_sig_mon,
                    j_sig_mask, new_pair_sig_mask))

                is_rewr = comp_sigratio(basis, j, new_basis_idx)
                is_rewr && break
            end
            if (j_sig_idx == basis_sig_idx &&
                div(monomial(basis.sigs[j]), basis_pair_sig_mon,
                    j_sig_mask, basis_pair_sig_mask))
                
                is_rewr = comp_sigratio(basis, j, i)
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
            if index(basis.sigs[j]) < basis_sig_idx
                is_rewr = div(leading_monomial(basis, basis_ht, j),
                              basis_pair_sig_mon,
                              basis.lm_masks[j], basis_pair_sig_mask)
                is_rewr && break
            end
        end
        is_rewr && continue
        
        # TODO: do we need to distinguish between top and bottom sig
        # TODO: this is wrong
        pair_deg = new_pair_sig_mon.deg + basis.degs[new_sig_idx]
        new_pair = if basis_sig_idx < new_sig_idx || (basis_sig_idx == new_sig_idx && comp_sigratio(basis, new_basis_idx, i))
                SPair(Sig(new_sigidx, new_pair_sig_mon),
                      Sig(basis_sig_idx, basis_pair_sig_mon),
                      new_pair_sig_mask, basis_pair_sig_mask,
                      new_basis_idx, i, pair_deg)
            else
                SPair(Sig(basis_sig_idx, basis_pair_sig_mon),
                      Sig(new_sigidx, new_pair_sig_mon),
                      basis_pair_sig_mask, new_pair_sig_mask,
                      i, new_basis_idx, pair_deg)
        end
        
        pairset.pairs[pairset.load + 1] = new_pair
        pairset.load += 1
    end
end
