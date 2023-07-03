# updating the pairset and basis

# add new reduced rows to basis/syzygies
# TODO: make a new symbol ht everytime?
function update_basis!(basis::Basis,
                       matrix::MacaulayMatrix,
                       pairset::Pairset{N},
                       symbol_ht::MonomialHashtable,
                       basis_ht::MonomialHashtable) where N

    add_indices = matrix.toadd
    # for now make sure that we add elements in signature order
    # TODO: shouldn't these be already sorted?
    sort!(add_indices,
          lt = (i1, i2) -> matrix.sig_order[i1] < matrix.sig_order[i2]) 

    @inbounds for i in add_indices
        # determine if row is zero
        row = matrix.rows[i]
        new_sig = matrix.sigs[i]
        new_idx = index(new_sig)
        new_sig_mon = monomial(new_sig)
        new_sig_mask = (new_idx, divmask(new_sig_mon, basis_ht.divmap,
                                         basis_ht.ndivbits))
        if isempty(row)
            # make sure we have enough space
            if basis.syz_load == basis.syz_size
                basis.syz_size *= 2
                resize!(basis.syz_sigs, basis.syz_size)
                resize!(basis.syz_masks, basis.syz_size)
            end
            
            # add new syz sig
            l = basis.syz_load + 1
            basis.syz_sigs[l] = new_sig
            basis.syz_masks[l] = new_sig_mask
            basis.syz_load += 1

            # kill pairs with known syz signature
            @inbounds for j in 1:pairset.load
                p = pairset.elems[j]
                index(p.top_sig) != new_idx && continue
                if divch(new_sig_mon, p.top_sig, new_sig_mask, p.top_sig_mask)
                    pairset.elems[j].top_index = 0
                end
                index(p.bot_sig) != new_idx && continue
                if divch(new_sig_mon, p.bot_sig, new_sig_mask, p.bot_sig_mask)
                    pairset.elems[i].top_index = 0
                end
            end

            # remove pairs that became rewriteable in previous loop
            remove_red_pairs!(pairset)
        else
            # make sure we have enough space
            if basis.basis_load == basis.basis_size
                basis.basis_size *= 2
                resize!(basis.sigs, basis.basis_size)
                resize!(sigmasks, basis.basis_size)
                resize!(sigratios, basis.basis_size)
                resize!(lm_masks, basis.basis_size)
                resize!(monomials, basis.basis_size)
                resize!(coefficients, basis.basis_size)
                resize!(is_red, basis.basis_size)
            end
            
            # add to basis hashtable
            insert_in_basis_hash_table_pivots!(row, basis_ht, symbol_ht)
            # TODO: check that this is correct
            lm = basis_ht.exponents[first(row)]

            # add everything to basis
            l = basis.basis_load + 1
            basis.sigs[l] = new_sig
            basis.sigmasks[l] = new_sig_mask
            basis.sigratios[l] = divide(new_sig_mon, lm)
            basis.lm_masks[l] = divmask(lm, basis_ht.divmap, basis_ht.ndivbits)
            basis.monomials[l] = row
            basis.coefficients[l] = matrix.coeffs[i]
            basis.basis_load = l

            # build new pairs
            update_pairset!(pairset, basis, basis_ht, l)
        end
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
    bmask = basis.sigmasks[new_basis_idx]
    @inbounds for i in 1:pairset.load
        p = pairset.elems[i]
        (iszero(p.top_index) || index(p.top_sig) != new_sig_idx) && continue
        if divch(new_sig_mon, monomial(p.top_sig), mask(bmask), p.top_sig_mask)
            if comp_sigratio(basis, new_basis_idx, p.top_index)
                pairset.elems[i].top_index = 0
                continue
            end
        end
        new_sig_idx != index(p.bot_sig) && continue
        if divch(new_sig_mon, monomial(p.bot_sig), mask(bmask), p.bot_sig_mask)
            if comp_sigratio(basis, new_basis_idx, p.bot_index)
                pairset.elems[i].top_index = 0
                continue
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
        all(iszero, mult_new_elem.exps) && error("no")
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
                is_rewr = divch(basis.syz_sigs[j], new_pair_sig_mon,
                              basis.syz_masks[j][2], new_pair_sig_mask) 
                is_rewr && break
            end
            if index(basis.syz_masks[j]) == basis_sig_idx
                is_rewr = divch(basis.syz_sigs[j], basis_pair_sig_mon,
                              basis.syz_masks[j][2], basis_pair_sig_mask) 
                is_rewr && break
            end
        end
        is_rewr && continue

        # check both pair sigs against basis sigs
        @inbounds for j in basis.basis_offset:basis.basis_load-1
            j == new_basis_idx && continue
            j_sig_idx = index(basis.sigmasks[j])
            (j_sig_idx != new_sig_idx && j_sig_idx != basis_sig_idx) && continue
            j_sig_mask = basis.sigmasks[j][2]
            if (j_sig_idx == new_sig_idx && divch(monomial(basis.sigs[j]),
                                                  new_pair_sig_mon,
                                                  j_sig_mask,
                                                  new_pair_sig_mask))
                is_rewr = comp_sigratio(basis, j, new_basis_idx)
                is_rewr && break
            end
            if (j_sig_idx == basis_sig_idx && divch(monomial(basis.sigs[j]),
                                                    basis_pair_sig_mon,
                                                    j_sig_mask,
                                                    basis_pair_sig_mask))
                is_rewr = comp_sigratio(basis, j, i)
                is_rewr && break
            end
        end
        is_rewr && continue
            
        # check both pair signatures against koszul syzygies
        # TODO: should we store the indices with the lm masks
        @inbounds for j in basis.basis_offset:(new_basis_idx-1)
            if index(basis.sigs[j]) < new_basis_idx
                is_rewr = divch(leading_monomial(basis, basis_ht, j),
                              new_pair_sig_mon,
                              basis.lm_masks[j], new_pair_sig_mask)
                is_rewr && break
            end
            if index(basis.sigs[j]) < basis_sig_idx
                is_rewr = divch(leading_monomial(basis, basis_ht, j),
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
                SPair((new_sig_idx, new_pair_sig_mon),
                      (basis_sig_idx, basis_pair_sig_mon),
                      new_pair_sig_mask, basis_pair_sig_mask,
                      new_basis_idx, i, pair_deg)
            else
                SPair((basis_sig_idx, basis_pair_sig_mon),
                      (new_sig_idx, new_pair_sig_mon),
                      basis_pair_sig_mask, new_pair_sig_mask,
                      i, new_basis_idx, pair_deg)
        end
        pairset.elems[pairset.load + 1] = new_pair
        pairset.load += 1
    end
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
    j = 1
    @inbounds for i in 1:pairset.load
        iszero(pairset.elems[i].top_index) && continue
        pairset.elems[j] = pairset.elems[i]
        j += 1
    end
    pairset.load = j 
end
