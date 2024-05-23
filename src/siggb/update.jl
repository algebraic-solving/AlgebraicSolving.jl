# updating the pairset and basis

# add new reduced rows to basis/syzygies
function update_siggb!(basis::Basis,
                       matrix::MacaulayMatrix,
                       pairset::Pairset{N},
                       symbol_ht::MonomialHashtable,
                       basis_ht::MonomialHashtable,
                       ind_order::IndOrder,
                       tags::Tags,
                       tr::Tracer,
                       vchar::Val{Char},
                       syz_queue::Vector{SyzInfo},
                       mod_ord::Symbol=:DPOT) where {N, Char}

    new_basis_c = 0
    new_syz_c = 0

    toadd = matrix.toadd[1:matrix.toadd_length]
    added_unit = false

    cofac_ins_inds = Int[]

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
            process_syzygy!(basis, new_sig, new_sig_mask, pairset, tr, mod_ord)
            push!(syz_queue, (basis.syz_load, Dict{SigIndex, Bool}()))

            # for cofactor insertion in ndeg computation
            if gettag(tags, index(new_sig)) in [:ndeg, :sat]
                push!(cofac_ins_inds, i)
            end
        else
            new_basis_c += 1
            coeffs = matrix.coeffs[i]
            parent_ind = matrix.parent_inds[i]
            added_unit = add_basis_elem!(basis, pairset, basis_ht, symbol_ht,
                                         row, coeffs,
                                         new_sig, new_sig_mask, parent_ind,
                                         tr, ind_order, tags, mod_ord)
        end
    end

    insert_syz_cofacs!(basis, basis_ht, matrix.sigs[cofac_ins_inds],
                       pairset, vchar, tr, tags, ind_order)

    if new_basis_c != 0 || new_syz_c != 0
        @info "$(new_basis_c) new, $(new_syz_c) zero"
    end

    return added_unit
end

function add_basis_elem!(basis::Basis{N},
                         pairset::Pairset,
                         basis_ht::MonomialHashtable,
                         symbol_ht::MonomialHashtable,
                         row::Vector{MonIdx},
                         coeffs::Vector{Coeff},
                         new_sig::Sig,
                         new_sig_mask::MaskSig,
                         parent_ind::Int,
                         tr::Tracer,
                         ind_order::IndOrder,
                         tags::Tags,
                         mod_ord::Symbol) where N


    # make sure we have enough space
    resize_basis!(basis)
    
    # add to basis hashtable
    insert_in_basis_hash_table_pivots!(row, basis_ht, symbol_ht)
    lm = basis_ht.exponents[first(row)]
    @debug "new lm $(lm)"
    lm_mask = divmask(lm, basis_ht.divmap, basis_ht.ndivbits)
    s = new_sig

    # check if we're adding a unit
    if length(row) == 1 && all(iszero, lm.exps)
        return true
    end

    # add everything to basis
    l = basis.basis_load + 1
    basis.sigs[l] = new_sig
    basis.sigmasks[l] = new_sig_mask
    new_sig_ratio = divide(lm, monomial(new_sig))
    basis.sigratios[l] = new_sig_ratio

    basis.lm_masks[l] = lm_mask
    basis.monomials[l] = row
    basis.coefficients[l] = coeffs

    basis.mod_rep_known[l] = falses(ind_order.max_ind)
    basis.mod_reps[l] = Vector{Polynomial}(undef, ind_order.max_ind)

    basis.is_red[l] = false

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

    # if an existing sig further reduced we dont need the old element
    if basis.sigs[parent_ind] == new_sig && parent_ind >= basis.basis_offset
        basis.is_red[parent_ind] = true
    end 

    basis.rewrite_nodes[l+1] = [-1, parent_ind+1]
    basis.basis_load = l

    # update tracer info
    store_basis_elem!(tr, new_sig, l, basis.basis_size)
    
    # build new pairs
    update_pairset!(pairset, basis, basis_ht, l, ind_order, tags, mod_ord)

    return false
end

function process_syzygy!(basis::Basis,
                         new_sig::Sig,
                         new_sig_mask::MaskSig,
                         pairset::Pairset,
                         tr::Tracer,
                         mod_ord::Symbol) 

    new_idx = index(new_sig_mask)
    
    new_sig_mon = monomial(new_sig)

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

    # add info to tracer
    store_syz!(tr)

    # kill pairs with known syz signature
    @inbounds for j in 1:pairset.load
        p = pairset.elems[j]
        cond = index(p.top_sig) == new_idx
        if cond && divch(new_sig_mon, monomial(p.top_sig),
                         new_sig_mask[2], p.top_sig_mask)
            pairset.elems[j].top_index = 0
        end
        cond = index(p.bot_sig) == new_idx && (mod_ord == :DPOT || index(p.bot_sig) == index(p.top_sig))
        if cond && divch(new_sig_mon, monomial(p.bot_sig),
                         new_sig_mask[2], p.bot_sig_mask)
            pairset.elems[j].top_index = 0
        end
    end

    # remove pairs that became rewriteable in previous loop
    remove_red_pairs!(pairset)
end

function insert_syz_cofacs!(basis::Basis{N},
                            basis_ht::MonomialHashtable{N},
                            syz_sigs::Vector{Sig{N}},
                            pairset::Pairset{N},
                            vchar::Val{Char},
                            tr::Tracer,
                            tags::Tags,
                            ind_order::IndOrder) where {N, Char}

    isempty(syz_sigs) && return
    cofacs = Polynomial[]
    mat_ind = length(tr.mats)
    syz_ind = index(first(syz_sigs))

    # we only use this in nondeg computation
    @assert length(unique([index(s) for s in syz_sigs])) == 1

    # construct cofactors of zero reduction and ins in hashtable
    for new_sig in syz_sigs
        new_idx = index(new_sig)
        @info "constructing module"
        cofac = construct_module(new_sig, basis,
                                 basis_ht,
                                 mat_ind, tr,
                                 vchar,
                                 ind_order,
                                 new_idx)
        if isempty(cofacs)
            push!(cofacs, cofac)
        else
            cofacs[1] = add_pols(cofac..., cofacs[1]...,
                                 vchar, rand(one(Coeff):Coeff(Char)))
            sort_poly!(cofac, by = midx -> basis_ht.exponents[midx],
                       lt = lt_drl, rev = true)
            normalize_cfs!(cofac[1], vchar)
            push!(cofacs, cofac)
        end
    end

    # insert cofactors in system
    tag = gettag(tags, index(first(syz_sigs)))
    new_f_idx = zero(SigIndex)
    for (i, cofac) in enumerate(cofacs)
        ord_ind = if tag == :sat
            sat_inds = findall(tag -> tag == :sat, tags)
            findmin(sat_ind -> ind_order.ord[sat_ind], sat_inds)[1]
        else
            ind_order.ord[syz_ind]
        end
        if isone(i)
            sort_poly!(cofac, by = midx -> basis_ht.exponents[midx],
                       lt = lt_drl, rev = true)
            normalize_cfs!(cofac[1], vchar)
        end
        new_tg = if tag == :sat
            isone(i) ? :fsatins : :satins
        else
            isone(i) ? :fndegins : :ndegins
        end
        new_ind = add_new_sequence_element!(basis, basis_ht, tr, cofac...,
                                            ind_order, ord_ind, pairset,
                                            tags,
                                            new_tg = new_tg)
    end
end

# construct all pairs with basis element at new_basis_idx
# and perform corresponding rewrite checks
# assumes DPOT
function update_pairset!(pairset::Pairset{N},
                         basis::Basis,
                         basis_ht::MonomialHashtable,
                         new_basis_idx::Int,
                         ind_order::IndOrder,
                         tags::Tags,
                         mod_ord::Symbol=:DPOT) where N

    new_sig_mon = monomial(basis.sigs[new_basis_idx])
    new_sig_idx = index(basis.sigs[new_basis_idx])

    # check existing pairs for rewriteability against element
    # at new_basis_idx and for koszul rewriteability
    @inbounds bmask = basis.sigmasks[new_basis_idx]
    @inbounds new_lm = leading_monomial(basis, basis_ht, new_basis_idx)
    @inbounds new_lm_msk = basis.lm_masks[new_basis_idx]
    @inbounds for i in 1:pairset.load
        p = pairset.elems[i]
        if (cmp_ind_str(new_sig_idx, index(p.top_sig), ind_order)
            && !are_incompat(new_sig_idx, index(p.top_sig), ind_order))
            if divch(new_lm, monomial(p.top_sig), new_lm_msk, p.top_sig_mask)
                pairset.elems[i].top_index = 0
                continue
            end
        end
        if !iszero(p.bot_index) && (mod_ord == :DPOT || index(p.bot_sig) != index(p.top_sig))
            if (cmp_ind_str(new_sig_idx, index(p.bot_sig), ind_order)
                && !are_incompat(new_sig_idx, index(p.top_sig), ind_order))
                if divch(new_lm, monomial(p.bot_sig), new_lm_msk, p.bot_sig_mask)
                    pairset.elems[i].top_index = 0
                    continue
                end
            end
        end
    end

    # kill pairs that became rewriteable in the previous round
    remove_red_pairs!(pairset)

    # resize pairset if needed
    num_new_pairs = new_basis_idx - 1
    resize_pairset!(pairset, num_new_pairs)

    new_sig_ratio = basis.sigratios[new_basis_idx]
    # pair construction loop
    @inbounds for i in basis.basis_offset:(new_basis_idx - 1)

        # ignore if redundant
        basis.is_red[i] && continue

        basis_sig_idx = index(basis.sigs[i])
        # dont build incompatible pairs 
        are_incompat(new_sig_idx, basis_sig_idx, ind_order) && continue

        basis_lm = leading_monomial(basis, basis_ht, i)
        mult_new_elem = lcm_div(new_lm, basis_lm)
        new_pair_sig_mon = mul(mult_new_elem, new_sig_mon)
        
        mult_basis_elem = lcm_div(basis_lm, new_lm)
        basis_pair_sig_mon = mul(mult_basis_elem, monomial(basis.sigs[i]))

        new_pair_sig = (new_sig_idx, new_pair_sig_mon)
        basis_pair_sig = (basis_sig_idx, basis_pair_sig_mon)

        new_sig_is_smaller = lt_pot(new_pair_sig, basis_pair_sig, ind_order)

        top_idx = new_sig_is_smaller ? basis_sig_idx : new_sig_idx
        
        # check if S-pair is singular
        new_pair_sig == basis_pair_sig && continue

        new_pair_sig_mask = divmask(new_pair_sig_mon,
                                    basis_ht.divmap,
                                    basis_ht.ndivbits)
        basis_pair_sig_mask = divmask(basis_pair_sig_mon,
                                      basis_ht.divmap,
                                      basis_ht.ndivbits)

        rewriteable_syz(basis, new_pair_sig,
                        new_pair_sig_mask, tags,
                        mod_ord == :DPOT || cmp_ind(basis_sig_idx, new_sig_idx, ind_order)) && continue
        rewriteable_syz(basis, basis_pair_sig,
                        basis_pair_sig_mask, tags,
                        mod_ord == :DPOT || cmp_ind(new_sig_idx, basis_sig_idx, ind_order)) && continue
        rewriteable_koszul(basis, basis_ht, new_pair_sig,
                           new_pair_sig_mask, ind_order, tags,
                           mod_ord == :DPOT || cmp_ind(basis_sig_idx, new_sig_idx, ind_order)) && continue
        rewriteable_koszul(basis, basis_ht, basis_pair_sig,
                           basis_pair_sig_mask, ind_order, tags,
                           mod_ord == :DPOT || cmp_ind(new_sig_idx, basis_sig_idx, ind_order)) && continue

        top_sig, top_sig_mask, top_index,
        bot_sig, bot_sig_mask, bot_index = begin
	    if new_sig_is_smaller
                basis_pair_sig, basis_pair_sig_mask, i,
                new_pair_sig, new_pair_sig_mask, new_basis_idx
            else
                new_pair_sig, new_pair_sig_mask, new_basis_idx,
                basis_pair_sig, basis_pair_sig_mask, i
            end 
        end
        
        pair_deg = new_pair_sig_mon.deg + basis.degs[new_sig_idx]
        new_pair = SPair(top_sig, bot_sig,
                         top_sig_mask, bot_sig_mask,
                         top_index, bot_index,
                         pair_deg)
        pairset.elems[pairset.load + 1] = new_pair
        pairset.load += 1
    end
end

function minimize!(basis::Basis{N},
                   pairset::Pairset{N},
                   tr::Tracer,
                   basis_ht::MonomialHashtable{N},
                   idx_bound::SigIndex,
                   ind_order::IndOrder,
                   tags::Tags) where N

    if !iszero(idx_bound)
        @info "minimizing up to $(ind_order.ord[idx_bound])"
    else
        @info "minimizing"
    end
    el_killed = 0

    sz = 10000
    min_data = Vector{Tuple{Monomial{N},
                            DivMask,
                            Int}}(undef, sz)

    j = 1
    @inbounds for i in basis.basis_offset:basis.basis_load
        bsi = index(basis.sigs[i])
        !iszero(idx_bound) && cmp_ind_str(idx_bound, bsi, ind_order) && continue
        gettag(tags, bsi) == :sat && continue
        if j > sz
            sz *= 2
            resize!(min_data, sz)
        end
        min_data[j] = (leading_monomial(basis, basis_ht, i),
                       basis.lm_masks[i], i)
        j += 1
    end

    l = j-1
    sort!(view(min_data, 1:l),
          by = x -> x[1],
          lt = (m1, m2) -> m1.deg <= m2.deg) 

    to_del = Int[]
    @inbounds for i in 1:l 
        lm, lm_msk, bind = min_data[i]
        basis.is_red[bind] && continue
        for j in i+1:l
            lm2, lm_msk2, bind2 = min_data[j]
            basis.is_red[bind2] && continue
            if divch(lm, lm2, lm_msk, lm_msk2)
                el_killed += 1
                basis.is_red[bind2] = true
                push!(to_del, bind2)
            end
        end
    end

    sort!(to_del)
    @info "$(el_killed) elements killed"
    garbage_collect!(basis, pairset,  tr, to_del)
end
