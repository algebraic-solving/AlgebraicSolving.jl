function select_normal!(pairset::Pairset{N},
                        basis::Basis{N},
                        matrix::MacaulayMatrix,
                        ht::MonomialHashtable,
                        symbol_ht::MonomialHashtable,
                        ind_order::IndOrder,
                        tags::Tags) where N

    # number of selected pairs
    npairs = 0
    min_pair_ind = 0
    deg = Exp(-1) 
    compat_ind = zero(SigIndex)
    dont_sel = Int[]
    for i in 1:pairset.load
        if deg == -1
            deg = pairset.elems[i].deg
            npairs += 1
            continue
        end
        if pairset.elems[i].deg == deg
            npairs += 1
            top_s_idx = index(pairset.elems[i].top_sig)
            if iszero(compat_ind) && gettag(tags, top_s_idx) == :sat
                compat_ind = top_s_idx
            elseif are_incompat(top_s_idx, compat_ind, ind_order)
                push!(dont_sel, i)
            end
        else
            break
        end
    end

    @info "selected $(npairs) pairs, degree $(deg)"
    added_to_matrix = 0

    # allocate matrix
    reinitialize_matrix!(matrix, npairs)
    skip = falses(npairs)

    l = 1
    dl = length(dont_sel)
    @inbounds for i in 1:npairs
        if skip[i]
            continue
        end

        if l <= dl && i == dont_sel[l]
            l += 1
            continue
        end

        pair = pairset.elems[i]
        # for each unique pair signature
        curr_top_sig = pair.top_sig

        rewriteable_syz(basis, curr_top_sig, pair.top_sig_mask, tags) && continue
        rewriteable_koszul(basis, ht,
                           curr_top_sig, pair.top_sig_mask,
                           ind_order, tags) && continue

        rewr_ind = find_canonical_rewriter(basis, pair.top_sig,
                                           pair.top_sig_mask)

        pair_with_rewr_ind = 0
        for j in i:npairs
            pair2 = pairset.elems[j]
            if pair2.top_sig == curr_top_sig
                skip[j] = true
                if pair2.top_index == rewr_ind
                    pair_with_rewr_ind = j
                end
            end
        end

        if !iszero(pair_with_rewr_ind)
            # take pair with non-rewr top signature
            pair = pairset.elems[pair_with_rewr_ind]
            
            # add if pair is a unit vector
            add_cond = iszero(pair.bot_index)

            # check if we have a pivot already
            if !add_cond
                mult = divide(monomial(curr_top_sig),
                              monomial(basis.sigs[rewr_ind]))
                lm = mul(mult, leading_monomial(basis, ht, rewr_ind))
                l_idx = find_in_hash_table(symbol_ht, lm)
                if !iszero(l_idx)
                    if matrix.pivot_size >= l_idx && !iszero(matrix.pivots[l_idx])
                        add_cond = matrix.sigs[matrix.pivots[l_idx]] != curr_top_sig
                    end
                end
            end
                
            # check if there is a non-rewriteable reducer
            if !add_cond
                mult = divide(monomial(curr_top_sig),
                              monomial(basis.sigs[rewr_ind]))
                lm = mul(mult, leading_monomial(basis, ht, rewr_ind))
                red_sig_mon = similar(symbol_ht.buffer)
                @inbounds for j in 1:N
                    red_sig_mon[j] = zero(Exp)
                end
                reducer_sig = (zero(SigIndex), monomial(red_sig_mon))
                reducer_ind = zero(SigIndex)
                
                @inbounds for j in 1:npairs
                    pair2 = pairset.elems[j]
                    pair2.bot_sig == curr_top_sig && continue
                    iszero(pair2.bot_index) && continue
                    are_incompat(index(pair2.bot_sig), compat_ind, ind_order) && continue
                    if iszero(reducer_ind) || lt_pot(pair2.bot_sig, reducer_sig, ind_order)
                        !lt_pot(pair2.bot_sig, curr_top_sig, ind_order) && continue
                        new_red = false
                        if !iszero(pair2.bot_index)
                            rewriteable(basis, ht, pair2.bot_index, pair2.bot_sig,
                                        pair2.bot_sig_mask, ind_order, tags) && continue
                            ind = pair2.bot_index
                            mult = divide(monomial(pair2.bot_sig),
                                          monomial(basis.sigs[ind]))
                            lm2 = mul(mult, leading_monomial(basis, ht, ind))
                            new_red = lm2 == lm
                        end
                        if new_red
                            reducer_sig = pair2.bot_sig
                            reducer_ind = pair2.bot_index
                        end
                    end
                end

                add_cond = !iszero(reducer_ind)
                if add_cond
                    mult = divide(monomial(reducer_sig),
                                  monomial(basis.sigs[reducer_ind]))
                    lead_idx = write_to_matrix_row!(matrix, basis,
                                                    reducer_ind,
                                                    symbol_ht, ht, mult,
                                                    reducer_sig)
                    # set pivot
                    resize_pivots!(matrix, symbol_ht)
                    matrix.pivots[lead_idx] = matrix.nrows
                end
            end

            # add row to matrix if any of the add conds is true
            if add_cond
                mult = divide(monomial(curr_top_sig),
                              monomial(basis.sigs[rewr_ind]))
                l_idx = write_to_matrix_row!(matrix, basis, rewr_ind,
                                             symbol_ht, ht, mult,
                                             curr_top_sig)
                added_to_matrix += 1
                if iszero(pair.bot_index)
                    matrix.toadd[matrix.toadd_length+1] = matrix.nrows
                    matrix.toadd_length += 1
                end
            end
        end
    end
                
    if !iszero(added_to_matrix)
        @info "$(added_to_matrix) non-rewriteable critical signatures added to matrix"
    end

    # remove selected pairs from pairset
    l = 1
    @inbounds for i in 1:(pairset.load-npairs+dl)
        if l <= dl && i == dont_sel[j]
            l += 1
            continue
        end
        pairset.elems[i] = pairset.elems[i+npairs-l+1]
    end
    pairset.load -= npairs
    return deg, compat_ind
end

function symbolic_pp!(basis::Basis{N},
                      matrix::MacaulayMatrix,
                      ht::MonomialHashtable,
                      symbol_ht::MonomialHashtable,
                      ind_order::IndOrder,
                      tags::Tags,
                      compat_ind::SigIndex=zero(SigIndex)) where N

    i = one(MonIdx)
    mult = similar(ht.buffer)
    mult2 = similar(ht.buffer)
    red_sig_mon = similar(ht.buffer)

    resize_pivots!(matrix, symbol_ht)

    # iterate over monomials in symbolic ht
    @inbounds while i <= symbol_ht.load
        found_reducer = false
        # skip if reducer already exists
        if !iszero(matrix.pivots[i])
            i += one(MonIdx)
            continue
        end

        red_ind = 0
        @inbounds for j in 1:N
            red_sig_mon[j] = zero(Exp)
        end
        mul_red_sig = (zero(SigIndex), monomial(red_sig_mon))

        # realloc matrix if necessary
        if matrix.size == matrix.nrows
            matrix.size *= 2
            resize!(matrix.rows, matrix.size)
            resize!(matrix.sigs, matrix.size)
            resize!(matrix.parent_inds, matrix.size)
            resize!(matrix.coeffs, matrix.size)
            resize!(matrix.toadd, matrix.size)
        end

        exp = symbol_ht.exponents[i]
        divm = symbol_ht.hashdata[i].divmask
        
        j = basis.basis_offset 
        @label target
        # find element in basis which divmask divides divmask of monomial
        @inbounds while j <= basis.basis_load && !basis.is_red[j] && !divch(basis.lm_masks[j], divm)
            j += 1
        end

        if j <= basis.basis_load
            @inbounds red_exp = leading_monomial(basis, ht, j)

            # actual divisibility check
            if !(divch!(mult2, exp, red_exp))
                j += 1
                @goto target
            end

            cand_sig = basis.sigs[j]
            if !iszero(compat_ind) && are_incompat(index(cand_sig), compat_ind, ind_order)
                j += 1
                @goto target
            end
            
            mul_cand_sig = (index(cand_sig),
                            mul(monomial(mult2), monomial(cand_sig)))
            cand_sig_mask = divmask(monomial(mul_cand_sig), ht.divmap,
                                    ht.ndivbits)

            # check if new reducer sig is smaller than possible previous
            if !iszero(red_ind) && lt_pot(mul_red_sig, mul_cand_sig, ind_order)
                j += 1
                @goto target
            end

            # check if reducer is rewriteable
            if rewriteable(basis, ht, j, mul_cand_sig, cand_sig_mask,
                           ind_order, tags)
                j += 1
                @goto target
            end

            # set reducer
            red_ind = j
            mul_red_sig = mul_cand_sig
            j += 1
            @inbounds for k in 1:N
                mult[k] = mult2[k]
            end
            @goto target
        end

        # write to matrix
        if !iszero(red_ind)
            mm = monomial(SVector(mult))
            if iszero(compat_ind) && gettag(tags, index(mul_red_sig)) == :sat
                compat_ind = index(mul_red_sig)
            end
            @inbounds lead_idx = write_to_matrix_row!(matrix, basis,
                                                      red_ind, symbol_ht,
                                                      ht, mm,
                                                      mul_red_sig)
            resize_pivots!(matrix, symbol_ht)
            matrix.pivots[lead_idx] = matrix.nrows
        end
        i += one(MonIdx)
    end
end

function finalize_matrix!(matrix::MacaulayMatrix,
                          symbol_ht::MonomialHashtable,
                          ind_order::IndOrder)
    
    # store indices into hashtable in a sorted way
    ncols = symbol_ht.load
    matrix.ncols = ncols

    col2hash = Vector{ColIdx}(undef, ncols)
    @inbounds for i in 1:ncols
        col2hash[i] = i
    end
    exps = symbol_ht.exponents
    function cmp(h1, h2)
        @inbounds e1 = exps[h1]
        @inbounds e2 = exps[h2]
        return !lt_drl(e1, e2)
    end
    sort!(col2hash, lt = cmp)
    matrix.col2hash = col2hash
    nc = matrix.ncols
    @inbounds matrix.pivots[1:nc] = matrix.pivots[1:nc][col2hash]

    if !iszero(matrix.nrows)
        @info "matrix of size $((matrix.nrows, matrix.ncols)), density $(@sprintf "%.2f" 100*sum((length).(matrix.rows[1:matrix.nrows]))/(matrix.nrows * matrix.ncols))%"
    end
    matrix.sig_order = Vector{Int}(undef, matrix.nrows)
    # sort signatures
    sortperm!(matrix.sig_order, matrix.sigs[1:matrix.nrows],
              lt = (sig1, sig2) -> lt_pot(sig1, sig2, ind_order))
end

# TODO: later to optimize: mem allocations for matrix
# helper functions
function initialize_matrix(::Val{N}) where {N}
    rows = Vector{Vector{MonIdx}}(undef, 0)
    pivots = Vector{Int}(undef, 0)
    pivot_size = 0
    sigs = Vector{Sig{N}}(undef, 0)
    parent_inds = Vector{Int}(undef, 0)
    sig_order = Vector{Int}(undef, 0)
    col2hash = Vector{ColIdx}(undef, 0)
    coeffs = Vector{Vector{Coeff}}(undef, 0)
    toadd = Vector{Int}(undef, 0)

    size = 0
    nrows = 0
    ncols = 0
    toadd_length = 0

    return MacaulayMatrix(rows, pivots, pivot_size,
                          sigs, parent_inds, sig_order,
                          col2hash, coeffs, size, nrows,
                          ncols, toadd, toadd_length)
end
    
# Refresh and initialize matrix for `npairs` elements
function reinitialize_matrix!(matrix::MacaulayMatrix, npairs::Int)
    matrix.size = 2 * npairs
    matrix.pivot_size = 2 * npairs
    resize!(matrix.rows, matrix.size)
    resize!(matrix.pivots, matrix.pivot_size)
    for i in 1:matrix.pivot_size
        matrix.pivots[i] = 0
    end
    resize!(matrix.sigs, matrix.size)
    resize!(matrix.parent_inds, matrix.size)
    resize!(matrix.coeffs, matrix.size)
    resize!(matrix.toadd, matrix.size)
    for i in 1:npairs
        matrix.toadd[i] = 0
    end
    return matrix
end

# resize pivots array if needed
@inline function resize_pivots!(matrix::MacaulayMatrix,
                                symbol_ht::MonomialHashtable)
    if matrix.pivot_size < symbol_ht.load 
        pv_size = matrix.pivot_size
        new_pv_size = 2 * (symbol_ht.load)
        resize!(matrix.pivots, new_pv_size)
        @inbounds for j in pv_size+1:new_pv_size 
            matrix.pivots[j] = 0
        end
        matrix.pivot_size = new_pv_size
    end
end
    
# helper function to write row to matrix
function write_to_matrix_row!(matrix::MacaulayMatrix,
                              basis::Basis,
                              basis_idx::Int,
                              symbol_ht::MonomialHashtable,
                              ht::MonomialHashtable,
                              mult::Monomial,
                              sig::Sig,
                              row_ind=matrix.nrows+1)

    hsh = Base.hash(mult)
    poly = basis.monomials[basis_idx]
    row = similar(basis.monomials[basis_idx])
    check_enlarge_hashtable!(symbol_ht, length(basis.monomials[basis_idx]))

    s = basis.sigs[basis_idx]
    lm = mul(mult, leading_monomial(basis, ht, basis_idx))

    @inbounds matrix.rows[row_ind] =
        insert_multiplied_poly_in_hash_table!(row, hsh, mult, poly,
                                              ht, symbol_ht)
    @inbounds matrix.coeffs[row_ind] = basis.coefficients[basis_idx]
    @inbounds matrix.sigs[row_ind] = sig
    @inbounds matrix.parent_inds[row_ind] = basis_idx
    if row_ind == matrix.nrows + 1
        matrix.nrows += 1
    end
    return first(matrix.rows[row_ind])
end
