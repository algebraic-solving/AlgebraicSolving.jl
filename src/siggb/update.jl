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
                       vchar::Val{Char}) where {N, Char}

    new_basis_c = 0
    new_syz_c = 0

    toadd = matrix.toadd[1:matrix.toadd_length]

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

            process_syzygy!(basis, basis_ht, pairset, new_sig, new_sig_mask,
                            tags, ind_order, tr, vchar)
        else
            new_basis_c += 1
            coeffs = matrix.coeffs[i]
            parent_ind = matrix.parent_inds[i]
            add_basis_elem!(basis, pairset, basis_ht, symbol_ht,
                            row, coeffs,
                            new_sig, new_sig_mask, parent_ind,
                            tr, ind_order, tags)
        end
    end
    if new_basis_c != 0 || new_syz_c != 0
        @info "$(new_basis_c) new, $(new_syz_c) zero"
    end
end

# add new reduced rows to basis/syzygies
function update_sigdecomp!(queue::Vector{Tuple{Basis{N}, Pairset{N}}},
                           basis::Basis,
                           matrix::MacaulayMatrix,
                           pairset::Pairset{N},
                           symbol_ht::MonomialHashtable,
                           basis_ht::MonomialHashtable,
                           ind_order::IndOrder,
                           tags::Tags,
                           tr::Tracer,
                           vchar::Val{Char}) where {N, Char}

    new_basis_c = 0
    new_syz_c = 0

    toadd = matrix.toadd[1:matrix.toadd_length]

    split_row_index = 0
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

            if !iszero(split_row_index) && gettag(tags, index(new_sig)) == :seq
                split_row_index = i
            else
                process_syzygy!(basis, basis_ht, pairset, new_sig, new_sig_mask,
                                tags, ind_order, tr, vchar)
            end
        else
            new_basis_c += 1
            coeffs = matrix.coeffs[i]
            parent_ind = matrix.parent_inds[i]
            add_basis_elem!(basis, pairset, basis_ht, symbol_ht,
                            row, coeffs,
                            new_sig, new_sig_mask, parent_ind,
                            tr, ind_order, tags)
        end
    end

    @inbounds if !iszero(split_row_index)
        # go through the splitting
        zd_sig = matrix.sigs[split_row_index]
        zd_ind = index(zd_sig)
        
        # initialize and fill new basis with
        basis_new = new_basis(basis.basis_size, basis.syz_size,
                              basis.input_size)

        # add every input element except for
        # the one in index zd_ind
        j = 1
        i = 1
        old_to_new_index = Vector{Int}(undef, basis.basis_load+1)
        incr_ind = ind -> ind == basis.input_load ? basis.basis_offset : ind+1
        while i < basis.basis_load
            s_ind = index(sig)

            # skip if index is the same as the zd
            if s_ind == zd_ind
                i = incr_ind(i)
                continue
            end

            # skip if index is > index of zd and not input el
            if ind_order.ord[s_ind] > ind_order.ord[zd_ind] && i > basis.input_load
                i = incr_ind(i)
                continue
            end
            
            # write into new basis
            overwrite_basis!(basis, basis_new, i, j)
            i <= basis.input_load ? basis_new.input_load += 1 : basis_new.basis_load += 1
            old_to_new_index[i+1] = j+1

            # fix rewrite tree
            if i <= basis.input_load
                push!(new_basis.rewrite_nodes[1], j)
                new_basis.rewrite_nodes[j+1] = [-1, 1]
            else
                rnodes = basis_new.rewrite_nodes[j+1]

                old_par_ind = rnodes[2]
                new_par_ind = old_to_new_index[old_parent_ind]

                # fix parent index
                rnodes[2] = new_par_ind
                
                # fix child index in parent
                par_nodes = basis_new.rewrite_nodes[new_par_ind]
                k = findfirst(ch -> ch == i+1, par_nodes[3:end])
                par_nodes[k+2] = j+1 
            end

            i = incr_ind(i)
        end

        # TODO: tbc...
    end
    
    if new_basis_c != 0 || new_syz_c != 0
        @info "$(new_basis_c) new, $(new_syz_c) zero"
    end
end

function add_basis_elem!(basis::Basis,
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
                         tags::Tags)

    # make sure we have enough space
    if basis.basis_load == basis.basis_size
        resize_basis!(basis)
    end
    
    # add to basis hashtable
    insert_in_basis_hash_table_pivots!(row, basis_ht, symbol_ht)
    lm = basis_ht.exponents[first(row)]
    s = new_sig

    # add everything to basis
    l = basis.basis_load + 1
    basis.sigs[l] = new_sig
    basis.sigmasks[l] = new_sig_mask
    new_sig_ratio = divide(lm, monomial(new_sig))
    basis.sigratios[l] = new_sig_ratio

    basis.lm_masks[l] = divmask(lm, basis_ht.divmap,
                                basis_ht.ndivbits)
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
        cofac_coeffs, cofac_mons = construct_module(new_sig, basis, mat_ind, tr, vchar,
                                                    Dict{Sig, Vector{Polynomial{N}}}(),
                                                    ind_order.max_ind, ind_order, new_idx)[new_idx]
        cofac_mons_hashed = [insert_in_hash_table!(basis_ht, mon) for mon in cofac_mons]

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
        ind = ind_order.max_ind + one(SigIndex)
        ind_order.max_ind = ind

        # update index order
        col_inds = findall(tag -> tag == :col, tags)
        filter!(col_ind -> ind_order.ord[col_ind] > ind_order.ord[new_idx], col_inds)
        if !isempty(col_inds)
            min_larger_ind, _ = findmin(col_ind -> ind_order.ord[col_ind], col_inds)
            @inbounds for i in eachindex(ind_order.ord)
                ord_i = ind_order.ord[i]
                if ord_i >= min_larger_ind
                    ind_order.ord[i] += one(SigIndex)
                end
            end
            push!(ind_order.ord, min_larger_ind)
        else
            push!(ind_order.ord, ind)
        end

        # update tags
        tags[ind] = :colins

        # add cofactor to input sequence
        if basis.input_load == basis.input_size
            make_room_new_input_el!(basis, tr)
        end
        lm = first(cofac_mons)
        lm_divm = divmask(lm, basis_ht.divmap, basis_ht.ndivbits)
        add_input_element!(basis, pairset, ind, cofac_mons_hashed,
                           cofac_coeffs, lm_divm, lm)
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
        dont_build_pair(new_sig_idx, basis_sig_idx, tags, ind_order) && continue
        dont_build_pair(basis_sig_idx, new_sig_idx, tags, ind_order) && continue

        basis_lm = leading_monomial(basis, basis_ht, i)
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
	    if lt_pot(basis_pair_sig, new_pair_sig, ind_order)
                new_pair_sig, new_pair_sig_mask, new_basis_idx,
                basis_pair_sig, basis_pair_sig_mask, i
            else
                basis_pair_sig, basis_pair_sig_mask, i,
                new_pair_sig, new_pair_sig_mask, new_basis_idx
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

function gettag(tags::Tags, i::SigIndex)
    return get(tags, i, :seq)
end

function dont_build_pair(ind1::SigIndex, ind2::SigIndex,
                         tags::Tags, ind_order::IndOrder)

    tag1 = gettag(tags, ind1)
    tag2 = gettag(tags, ind2)

    tag1 == tag2 == :col && return true

    if tag1 == :colins
        ind1 == ind2 && return true
        if cmp_ind(ind2, ind1, ind_order)
            return true
        end
    end
    return false
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

function resize_basis!(basis::Basis)
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

function make_room_new_input_el!(basis::Basis,
                                 tr::Tracer)

    # this whole block just shifts the basis to the right
    # to make room for new input elements
    @inbounds begin
        basis.input_size *= 2
        
        shift = basis.input_size
        old_offset = basis.basis_offset
        basis.basis_offset += shift
        basis.basis_load += shift
        if basis.basis_load >= basis.basis_size
            resize_basis!(basis)
        end

        if tr.load + shift >= tr.size
            tr.size *= 2
            resize!(tr.mats, tr.size)
            resize!(tr.basis_ind_to_mat, tr.size)
        end

        # adjusts rewrite nodes at the start
        for i in 1:basis.input_load
            basis.rewrite_nodes[i+1][3] += shift
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
            for j in 3:rnodes[1]
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
    end
end
