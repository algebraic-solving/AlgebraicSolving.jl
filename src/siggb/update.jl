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
                       syz_queue::Vector{Int}) where {N, Char}

    new_basis_c = 0
    new_syz_c = 0

    toadd = matrix.toadd[1:matrix.toadd_length]
    added_unit = false

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

            # check if syzygy could be interesting to split with
            if all(iszero, new_sig_mon.exps)
                if gettag(tags, new_idx) == :col
                    return true
                end
                # mark input element as redundant
                basis.is_red[new_idx] = true
            elseif gettag(tags, new_idx) == :split
                push!(syz_queue, basis.syz_load+1)
            end

            process_syzygy!(basis, basis_ht, pairset, new_sig,
                            new_sig_mask,
                            tags, ind_order, tr, vchar)
        else
            new_basis_c += 1
            coeffs = matrix.coeffs[i]
            parent_ind = matrix.parent_inds[i]
            added_unit = add_basis_elem!(basis, pairset, basis_ht, symbol_ht,
                                         row, coeffs,
                                         new_sig, new_sig_mask, parent_ind,
                                         tr, ind_order, tags)
        end
    end

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
                         tags::Tags) where N


    # make sure we have enough space
    resize_basis!(basis)
    
    # add to basis hashtable
    insert_in_basis_hash_table_pivots!(row, basis_ht, symbol_ht)
    lm = basis_ht.exponents[first(row)]
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

    basis.is_red[l] = false

    if gettag(tags, index(new_sig)) != :colins
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
    end

    basis.rewrite_nodes[l+1] = [-1, parent_ind+1]
    basis.basis_load = l

    # update tracer info
    if tr.load >= tr.size
        tr.size *= 2
        resize!(tr.mats, tr.size)
        resize!(tr.basis_ind_to_mat, tr.size)
    end
    @inbounds tr.basis_ind_to_mat[l] = length(tr.mats)
    tr.load += 1

    # build new pairs
    update_pairset!(pairset, basis, basis_ht, l, ind_order, tags)

    return false
end

function process_syzygy!(basis::Basis{N},
                         basis_ht::MonomialHashtable,
                         pairset::Pairset,
                         new_sig::Sig,
                         new_sig_mask::MaskSig,
                         tags::Tags,
                         ind_order::IndOrder,
                         tr::Tracer,
                         vchar::Val{Char}) where {N, Char}

    new_idx = index(new_sig_mask)
    tag = gettag(tags, new_idx)
    
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
    push!(tr.syz_ind_to_mat, length(tr.mats)) 

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
    if tag == :col
        @info "inserting cofactor from colon ideal computation"
        # construct cofactor of zero reduction and ins in hashtable
        mat_ind = length(tr.mats)
        @info "constructing module"
        cofac_coeffs, cofac_mons = construct_module(new_sig, basis,
                                                    mat_ind, tr,
                                                    vchar,
                                                    ind_order.max_ind,
                                                    ind_order, new_idx)[new_idx]
        @info "done"
        cofac_mons_hashed = [insert_in_hash_table!(basis_ht, mon)
                             for mon in cofac_mons]

        # normalize coefficients
        inver = inv(first(cofac_coeffs), vchar)
        @inbounds for i in eachindex(cofac_coeffs)
            if isone(i)
                cofac_coeffs[i] = one(Coeff)
                continue
            end
            cofac_coeffs[i] = mul(inver, cofac_coeffs[i], vchar)
        end

        # sig index for cofactor
        ind = SigIndex(basis.input_load + 1)

        # update index order
        col_inds = findall(tag -> tag == :col, tags)
        filter!(col_ind -> ind_order.ord[col_ind] > ind_order.ord[new_idx],
                col_inds)
        if !isempty(col_inds)
            ord_ind, _ = findmin(col_ind -> ind_order.ord[col_ind], col_inds)
        else
            ord_ind = ind_order.max_ind + one(SigIndex)
        end
        ins_index!(ind_order, ord_ind)
        ind_order.incompat[(new_idx, ind)] = true

        # update tags
        tags[ind] = :colins

        # add cofactor to input sequence
        if basis.input_load == basis.input_size
            make_room_new_input_el!(basis, pairset, tr)
        end

        lm = first(cofac_mons)

        lm_divm = divmask(lm, basis_ht.divmap, basis_ht.ndivbits)
        add_input_element!(basis, ind, cofac_mons_hashed,
                           cofac_coeffs, lm_divm, lm)
        add_unit_pair!(basis, pairset, ind, lm.deg)
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
                         tags::Tags) where N


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

        # ignore if redundant
        basis.is_red[i] && continue

        basis_sig_idx = index(basis.sigs[i])
        # dont build some pairs if one of the elements is inserted
        # during colon ideal computation
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
        gettag(tags, top_idx) == :colins && continue
        
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
                        new_pair_sig_mask, tags) && continue
        rewriteable_syz(basis, basis_pair_sig,
                        basis_pair_sig_mask, tags) && continue

        # check both pair signatures against koszul syzygies
        rewriteable_koszul(basis, basis_ht, new_pair_sig,
                           new_pair_sig_mask, ind_order, tags) && continue
        rewriteable_koszul(basis, basis_ht, basis_pair_sig,
                           basis_pair_sig_mask, ind_order, tags) && continue

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

# helper functions
function leading_monomial(basis::Basis,
                          basis_ht::MonomialHashtable,
                          i)

    return basis_ht.exponents[first(basis.monomials[i])]
end

function gettag(tags::Tags, i::Integer)
    return get(tags, SigIndex(i), :seq)
end

function sort_pairset_by_degree!(pairset::Pairset, from::Int, sz::Int)
    ordr = Base.Sort.ord(isless, p -> p.deg, false, Base.Sort.Forward)
    sort!(pairset.elems, from, from+sz, def_sort_alg, ordr) 
end

# insert new index
function ins_index!(ind_order::IndOrder,
                    new_ord_ind::Integer)

    @inbounds for i in eachindex(ind_order.ord)
        if ind_order.ord[i] >= new_ord_ind
            ind_order.ord[i] += one(SigIndex)
        end
    end
    push!(ind_order.ord, SigIndex(new_ord_ind))
    ind_order.max_ind += 1
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

#---------------- resize functions --------------------#

function resize_basis!(basis::Basis)
    if basis.basis_load >= basis.basis_size
        basis.basis_size *= 2
        resize!(basis.sigs, basis.basis_size)
        resize!(basis.sigmasks, basis.basis_size)
        resize!(basis.sigratios, basis.basis_size)
        resize!(basis.rewrite_nodes, basis.basis_size)
        resize!(basis.lm_masks, basis.basis_size)
        resize!(basis.monomials, basis.basis_size)
        resize!(basis.coefficients, basis.basis_size)
        resize!(basis.is_red, basis.basis_size)
        resize!(basis.mod_rep_known, basis.basis_size)
        resize!(basis.mod_reps, basis.basis_size)
    end
end

function resize_syz!(basis::Basis)
    if basis.syz_load >= basis.syz_size
        basis.syz_size *= 2
        resize!(basis.syz_sigs, basis.syz_size)
        resize!(basis.syz_masks, basis.syz_size)
    end
end

function make_room_new_input_el!(basis::Basis,
                                 pairset::Pairset,
                                 tr::Tracer)

    # this whole block just shifts the basis to the right
    # to make room for new input elements
    @inbounds begin
        basis.input_size *= 2
        
        shift = basis.input_size
        old_offset = basis.basis_offset
        basis.basis_offset += shift
        basis.basis_load += shift
        resize_basis!(basis)

        if tr.load + shift >= tr.size
            tr.size *= 2
            resize!(tr.mats, tr.size)
            resize!(tr.basis_ind_to_mat, tr.size)
        end

        # adjusts rewrite nodes at the start
        for i in 1:basis.input_load
            if basis.rewrite_nodes[i+1][1] >= 0
                basis.rewrite_nodes[i+1][3] += shift
            end
        end

        for i in basis.basis_load:-1:basis.basis_offset
            basis.sigs[i] = basis.sigs[i-shift]
            basis.sigmasks[i] = basis.sigmasks[i-shift]
            basis.sigratios[i] = basis.sigratios[i-shift]
            basis.lm_masks[i] = basis.lm_masks[i-shift]
            basis.monomials[i] = basis.monomials[i-shift]
            basis.coefficients[i] = basis.coefficients[i-shift]
            basis.is_red[i] = basis.is_red[i-shift]

            # adjust tracer
            tr.basis_ind_to_mat[i] = tr.basis_ind_to_mat[i-shift]
            
            # adjust rewrite tree
            rnodes = basis.rewrite_nodes[i-shift+1]
            if rnodes[2] >= old_offset + 1
                rnodes[2] += shift
            end
            for j in 3:3+rnodes[1]
                rnodes[j] += shift
            end
            basis.rewrite_nodes[i+1] = rnodes
        end

        # adjust tracer
        for mat in tr.mats
            for i in keys(mat.rows)
                v = mat.rows[i]
                if v[2] >= old_offset
                    mat.rows[i] = (v[1], v[2] + shift)
                end
            end
        end

        # adjust pairset
        for i in 1:pairset.load
            p = pairset.elems[i]
            if p.top_index >= old_offset
                pairset.elems[i].top_index += shift
            end
            if p.bot_index >= old_offset
                pairset.elems[i].bot_index += shift
            end
        end
    end
end


# # this membership test is temporary
# does_div = false
# for j in basis.basis_offset:basis.basis_load
#     tg = gettag(tags, index(basis.sigs[j]))
#     tg == :col && continue
#     lm_msk = basis.lm_masks[j]
#     lm = basis_ht.exponents[first(basis.monomials[j])]
#     if divch(lm, new_sig_mon,
#              lm_msk, mask(new_sig_mask))
#         does_div = true
#         break
#     end
# end

# if !does_div
#     @info "found something to split with"
#     mat_ind = length(tr.mats)
#     cofac_coeffs, cofac_mons = construct_module(new_sig, basis, mat_ind, tr, vchar,
#                                                 ind_order.max_ind,
#                                                 ind_order, new_idx)[new_idx]
#     cofac_mons_hashed = [insert_in_hash_table!(basis_ht, mon) for mon in cofac_mons]
#     return true, cofac_coeffs, cofac_mons_hashed, new_idx
# else
#     @info "cofac contained in ideal, not splitting"
# end
