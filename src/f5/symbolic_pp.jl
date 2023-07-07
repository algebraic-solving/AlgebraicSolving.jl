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

    k = 1
    @inbounds for i in 1:npairs
        skip[i] && continue
        pair = pairset.elems[i]
        # for each unique pair signature
        curr_top_sig = pair.top_sig
        reducer_sig = pair.bot_sig
        reducer_ind = pair.bot_index

        # add row to be reduced to matrix
        mult = divide(monomial(pair.top_sig),
                      monomial(basis.sigs[pair.top_index]))
        write_to_matrix_row!(matrix, basis, pair.top_index, symbol_ht,
                             ht, mult, pair.top_sig) 

        # mark it to be added later
        matrix.toadd[k] = matrix.nrows
        k += 1

        # find the minimal top reducing bottom signature
        # input elements are stored as pairs with bot_index = 0
        if !iszero(reducer_ind)
            # add reducer row
            @inbounds for j in (i+1):npairs
                pair2 = pairset.elems[j]
                if pair2.top_sig == curr_top_sig
                    skip[j] = true
                    if lt_pot(pair2.bot_sig, reducer_sig)
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

            # set pivot
            resize_pivots!(matrix, symbol_ht)
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

    # iterate over monomials in symbolic ht
    @inbounds while i <= symbol_ht.load
        # skip if reducer already exists
        if !iszero(matrix.pivots[i])
            i += one(MonIdx)
            continue
        end

        # realloc matrix if necessary
        if matrix.size == matrix.nrows
            matrix.size *= 2
            resize!(matrix.rows, matrix.size)
            resize!(matrix.sigs, matrix.size)
            resize!(matrix.coeffs, matrix.size)
        end

        exp = symbol_ht.exponents[i]
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
            if !(divch!(mult, exp, red_exp))
                j += 1
                @goto target
            end

            # set reducer
            red_ind = j
            @inbounds red_sig = basis.sigs[j]

            # now that we found a reducer, we start looking for a better one
            @label target2
            j += 1
            @inbounds while j <= basis.basis_load && !divch(basis.lm_masks[j], divm)
                j += 1 
            end

            @inbounds if j <= basis.basis_load
                cand_sig = basis.sigs[j]
                cand_index = index(cand_sig)

                # skip if index is larger than current reducer
                if cand_index > index(red_sig)
                    @goto target2
                end

                # actual divisibility check
                div_flag = true
                @inbounds cand_exp = leading_monomial(basis, ht, j)
                if !(divch!(mult2, exp, cand_exp))
                    @goto target2
                end

                # check if new candidate reducer has smaller signature
                if (cand_index < index(red_sig) ||
                    lt_drl(mul(monomial(mult2), monomial(cand_sig)),
                           mul(monomial(mult), monomial(red_sig))))

                    mult = mult2
                    red_ind = j
                    red_sig = cand_sig
                    @goto target2
                end

                # check if new candidate rewrites reducer
                # TODO: in theory the following is correct?
                if (divch(monomial(cand_sig),
                          mul(monomial(mult), monomial(red_sig))) &&
                    comp_sigratio(basis, red_ind, j))
                    mult = mult2
                    red_ind = j
                    red_sig = cand_sig
                    @goto target2
                end
            end
            mm = monomial(SVector(mult))
            mul_red_sig = (index(red_sig), mul(mm, monomial(red_sig)))
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

    # set pivots correctly
    # println("PIVOTS: ")
    # println((Int).(matrix.pivots[1:matrix.ncols]))
    # @inbounds for i in 1:matrix.ncols
    #     matrix.pivots[hash2col[i]] = matrix.pivots[i]
    # end
    # @inbounds matrix.pivots = matrix.pivots[hash2col]
    # println("PIVOTS: ")
    # println((Int).(matrix.pivots[1:matrix.ncols]))

    # sort signatures
    @info "matrix of size $((matrix.nrows, matrix.ncols)), density $(sum((length).(matrix.rows[1:matrix.nrows]))/(matrix.nrows * matrix.ncols))"
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

    return MacaulayMatrix(rows, pivots, pivot_size,
                          sigs, sig_order, col2hash,
                          coeffs, size, nrows, ncols,
                          toadd)
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
    resize!(matrix.toadd, npairs)
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
                              sig::Sig)

    hsh = Base.hash(mult)
    row_ind = matrix.nrows + 1
    poly = basis.monomials[basis_idx]
    row = similar(basis.monomials[basis_idx])
    check_enlarge_hashtable!(symbol_ht, length(basis.monomials[basis_idx]))
    @inbounds matrix.rows[row_ind] =
        insert_multiplied_poly_in_hash_table!(row, hsh, mult, poly,
                                              ht, symbol_ht)
    @inbounds matrix.coeffs[row_ind] = basis.coefficients[basis_idx]
    @inbounds matrix.sigs[row_ind] = sig
    matrix.nrows += 1
    return first(matrix.rows[row_ind])
end
