function select_normal!(pairset::Pairset{N},
                        basis::Basis{N},
                        matrix::MacaulayMatrix,
                        ht::MonomialHashtable,
                        symbol_ht::MonomialHashtable) where N

    # sort pairset
    sort_pairset_by_degree!(pairset, 1, pairset.load-1)

    # number of selected pairs
    npairs = 0
    deg = zero(Exp) 
    for i in 1:pairset.load
        if iszero(deg)
            deg = pairset.elems[i].deg
            npairs += 1
            continue
        end
        if pairset.elems[i].deg == deg
            npairs += 1
        else
            break
        end
    end

    @info "selected $(npairs) pairs, degree $(deg)"

    # allocate matrix
    reinitialize_matrix!(matrix, npairs)
    skip = falses(npairs)

    @inbounds for i in 1:npairs
        if skip[i]
            continue
        end

        pair = pairset.elems[i]
        # for each unique pair signature
        curr_top_sig = pair.top_sig
        reducer_sig = pair.bot_sig
        reducer_ind = pair.bot_index

        # add row to be reduced to matrix
        mult = divide(monomial(pair.top_sig),
                      monomial(basis.sigs[pair.top_index]))
        
        l_idx = write_to_matrix_row!(matrix, basis, pair.top_index, symbol_ht,
                                     ht, mult, pair.top_sig) 
        sig = pair.top_sig
        lm = symbol_ht.exponents[l_idx]

        # mark it to be added later
        if iszero(reducer_ind)
            matrix.toadd[matrix.toadd_length+1] = matrix.nrows
            matrix.toadd_length += 1
        end

        @inbounds for j in (i+1):npairs
            pair2 = pairset.elems[j]
            if pair2.top_sig == curr_top_sig
                skip[j] = true
            end
        end

        # find the minimal top reducing bottom signature
        # input elements are stored as pairs with bot_index = 0
        resize_pivots!(matrix, symbol_ht)
        if !iszero(reducer_ind) && iszero(matrix.pivots[l_idx])

            @inbounds for j in (i+1):npairs
                pair2 = pairset.elems[j]
                if lt_pot(pair2.bot_sig, reducer_sig)
                    new_red = false
                    if pair2.top_sig == curr_top_sig
                        new_red = true
                    elseif !iszero(pair2.bot_index)
                        ind = pair2.bot_index
                        mult = divide(monomial(pair2.bot_sig),
                                      monomial(basis.sigs[ind]))
                        lm = mul(mult, leading_monomial(basis, ht, ind))
                        new_red = lm == symbol_ht.exponents[l_idx]
                    end
                    if new_red
                        reducer_sig = pair2.bot_sig
                        reducer_ind = pair2.bot_index
                    end
                end
            end

            mult = divide(monomial(reducer_sig),
                          monomial(basis.sigs[reducer_ind]))
            lead_idx = write_to_matrix_row!(matrix, basis, reducer_ind,
                                            symbol_ht, ht, mult,
                                            reducer_sig)
            sig = reducer_sig
            lm = symbol_ht.exponents[lead_idx]

            # set pivot
            matrix.pivots[lead_idx] = matrix.nrows
        end
    end

    # remove selected pairs from pairset
    @inbounds for i in 1:(pairset.load-npairs)
        pairset.elems[i] = pairset.elems[i+npairs]
    end
    pairset.load -= npairs
    resize_pivots!(matrix, symbol_ht)
end

function symbolic_pp!(basis::Basis{N},
                      matrix::MacaulayMatrix,
                      ht::MonomialHashtable,
                      symbol_ht::MonomialHashtable) where N

    i = one(MonIdx)
    mult = similar(ht.buffer)
    mult2 = similar(ht.buffer)
    red_sig_mon = similar(ht.buffer)

    # iterate over monomials in symbolic ht
    @inbounds while i <= symbol_ht.load
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
            resize!(matrix.coeffs, matrix.size)
            resize!(matrix.toadd, matrix.size)
        end

        exp = symbol_ht.exponents[i]
        # println("checking $(exp.exps)")
        divm = symbol_ht.hashdata[i].divmask
        
        j = basis.basis_offset 
        @label target
        # find element in basis which divmask divides divmask of monomial
        @inbounds while j <= basis.basis_load && !divch(basis.lm_masks[j], divm)
            j += 1 
        end

        if j <= basis.basis_load
            @inbounds red_exp = leading_monomial(basis, ht, j)

            # actual divisibility check
            if !(divch!(mult2, exp, red_exp))
                j += 1
                @goto target
            end

            # check if new reducer sig is smaller than possible previous
            cand_sig = basis.sigs[j]
            mul_cand_sig = (index(cand_sig),
                            mul(monomial(mult2), monomial(cand_sig)))
            if !iszero(red_ind) && lt_pot(mul_red_sig, mul_cand_sig)
                j += 1
                @goto target
            end

            # check if reducer is rewriteable
            cand_sig_mask = divmask(monomial(mul_cand_sig), ht.divmap,
                                    ht.ndivbits)
            if rewriteable(basis, ht, j, mul_cand_sig, cand_sig_mask)
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
                          symbol_ht::MonomialHashtable)
    
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

    # sort signatures
    @info "matrix of size $((matrix.nrows, matrix.ncols)), density $(@sprintf "%.2f" sum((length).(matrix.rows[1:matrix.nrows]))/(matrix.nrows * matrix.ncols))"
    matrix.sig_order = Vector{Int}(undef, matrix.nrows)
    sortperm!(matrix.sig_order, matrix.sigs[1:matrix.nrows],
              lt = (sig1, sig2) -> lt_pot(sig1, sig2))
end

# TODO: later to optimize: mem allocations for matrix
# helper functions
function initialize_matrix(::Val{N}) where {N}
    rows = Vector{Vector{MonIdx}}(undef, 0)
    pivots = Vector{Int}(undef, 0)
    pivot_size = 0
    sigs = Vector{Sig{N}}(undef, 0)
    sig_order = Vector{Int}(undef, 0)
    col2hash = Vector{ColIdx}(undef, 0)
    coeffs = Vector{Vector{Coeff}}(undef, 0)
    toadd = Vector{Int}(undef, 0)

    size = 0
    nrows = 0
    ncols = 0
    toadd_length = 0

    return MacaulayMatrix(rows, pivots, pivot_size,
                          sigs, sig_order, col2hash,
                          coeffs, size, nrows, ncols,
                          toadd, toadd_length)
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
    # println("row $((Int(index(sig)), monomial(sig).exps)), $(lm.exps)")

    @inbounds matrix.rows[row_ind] =
        insert_multiplied_poly_in_hash_table!(row, hsh, mult, poly,
                                              ht, symbol_ht)
    @inbounds matrix.coeffs[row_ind] = basis.coefficients[basis_idx]
    @inbounds matrix.sigs[row_ind] = sig
    if row_ind == matrix.nrows + 1
        matrix.nrows += 1
    end
    return first(matrix.rows[row_ind])
end
