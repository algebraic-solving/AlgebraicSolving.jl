# updating the pairset and basis

# for readability
function leading_monomial(basis::POTBasis,
                          basis_ht::MonomialHashtable,
                          i)

    return basis_ht.exponents[first(basis.monomials[i])]
end

# construct all pairs with basis element at new_basis_idx
# and perform corresponding rewrite checks
function update_pairset!(pairset::Pairset{POTSPair{N}},
                         basis::POTBasis,
                         basis_ht::MonomialHashtable,
                         new_basis_idx::Int) where N


    # resize pairset if needed
    pair_size = length(pairset.pairs)
    num_new_pairs = new_basis_idx - 1
    if pairset.load + num_new_pairs >= pair_size
        resize!(pairset.pairs, max(2 * pair_size, pair.load - num_new_pairs))
    end

    new_sig_mon = monomial(basis.sigs[new_basis_idx])

    # check existing pairs for rewriteability against element
    # at new_basis_idx
    bmask = basis.sigmasks[new_basis_idx]
    @inbounds for i in 1:(pairset.load)
        p = pairset.pairs[i]
        iszero(p.top_index) && continue
        if div(new_sig_mon, p.top_sig, bmask, p.top_sig_mask)
            pairset.pairs[i].top_index = 0
            continue
        end
        basis.curr_indx != index(p.bot_sig) && continue
        if div(new_sig_mon, p.bot_sig, bmask, p.bot_sig_mask)
            pairset.pairs[i].top_index = 0
            continue
        end
    end

    new_lm = basis_ht.exponents[basis.lms[new_basis_idx]]
    # pair construction loop
    @inbounds for i in 1:(new_basis_idx - 1)
        basis_lm = basis_ht.exponents[basis.lms[i]]
        ind = index(basis.sigs[i])

        mult_new_elem = lcm_div(new_lm, basis_lm)
        new_pair_sig_mon = mul(mult_new_elem, new_sig_mon)

        mult_basis_elem = lcm_div(basis_lm, new_lm)
        basis_pair_sig_mon = mul(mult_basis_elem, monomial(basis.sigs[i]))

        # check if S-pair is singular
        basis.curr_indx == ind && new_pair_sig_mon == basis_pair_sig_mon && continue

        new_pair_sig_mask = divmask(new_pair_sig_mon,
                                    basis_ht.divmap,
                                    basis_ht.ndivbits)
        basis_pair_sig_mask = divmask(basis_pair_sig_mon,
                                      basis_ht.divmap,
                                      basis_ht.ndivbits)

        # check both pair signatures against koszul syzygies
        exc_ind = false
        is_rewr = false
        for j in 1:(basis.curr_indx_start-1)
            is_rewr = div(leading_monomial(basis, basis_ht, j),
                          new_pair_sig_mon,
                          basis.lm_masks[j], new_pair_sig_mask)
            is_rewr && break
            if !exc_ind
                if j >= basis.index_cutoffs[ind+1]
                    exc_ind = true
                    continue
                end
                is_rewr = div(leading_monomial(basis, basis_ht, j),
                              basis_pair_sig_mon,
                              basis.lm_masks[j], basis_pair_sig_mask)
                is_rewr && break
            end
        end
        is_rewr && continue

        # check both pair sigs against non-trivial syzygies
        rewr_syz((basis.curr_indx, new_pair_sig_mon),
                 new_pair_sig_mask, basis) && continue
        rewr_syz((ind, basis_pair_sig_mon),
                 basis_pair_sig_mask, basis) && continue
        
        # check multiplied signature of basis element against basis sigs
        is_rewr = false
        ind_cuts = basis.index_cutoffs
        start_ind = ind_cuts[ind]
        end_ind = length(ind_cuts) == ind ? ind_cuts[ind+1]-1 : basis.basis_load
        @inbounds for j in start_ind:end_ind
            is_rewr = div(monomial(basis.sigs[j]), basis_pair_sig_mon,
                          basis.sigmasks[j], basis_pair_sig_mask)
            is_rewr && break
        end
        is_rewr && continue
        
        new_pair = if index(basis.sigs[i]) == basis.curr_indx && lt_drl(new_sig_mon, basis_pair_sig_mon)
            POTSPair(new_sig_mon, Sig(ind, basis_pair_sig_mon),
                     new_pair_sig_mask, basis_pair_sig_mask, new_basis_idx, i)
        else
            POTSPair(basis_pair_sig_mon, new_sig_mon,
                     basis_pair_mon_mask, new_pair_sig_mask, i, new_basis_idx)
        end
            
        pairset.pairs[pairset.load + 1] = new_pair
        pairset.load += 1
    end
end

@inline function rewr_syz(sig::Sig{N},
                          sigmask::DivMask,
                          basis::POTBasis{N}) where N

    sigmon = monomial(sig)
    ind = index(sig)
    syz_ind_cuts = basis.syz_index_cutoffs
    start_ind = syz_ind_cuts[ind]
    end_ind = length(syz_ind_cuts) == ind ? syz_ind_cuts[ind+1]-1 : basis.syz_load
    @inbounds for j in start_ind:end_ind
        div(basis.syz_sigs[j], sigmon,
            basis.syz_masks[j], sigmask) && return true
    end

    return false 
end
