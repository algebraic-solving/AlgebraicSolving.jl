# TODO: enforce some type restrictions
function select_normal!(pairset::Pairset{SPair{N}},
                        basis::Basis{N, C},
                        matrix::MacaulayMatrix{C},
                        ht::MonomialHashtable,
                        symbol_ht::MonomialHashtable) where {N, C}

    # number of selected pairs
    npairs = 0
    deg = zero(Exp) 
    for i in 1:pairset.load
        if iszero(deg)
            deg = pairset.pairs[i].deg
            npairs += 1
        end
        if pairset.pairs[i].deg == deg
            npairs += 1
        else
            break
        end
    end

    # allocate matrix
    reinitialize_matrix!(matrix, npairs)
    skip = falses(npairs)

    @inbounds for i in 1:npairs
        skip[i] && continue
        pair = pairset.pairs[i]
        # for each unique pair signature
        curr_top_sig = pair.top_sig
        reducer_sig = pair.bot_sig

        # find the minimal top reducing bottom signature
        @inbounds for j in (i+1):npairs
            pair2 = pairset.pairs[j]
            if pair2.top_sig == curr_top_sig
                skip[j] = true
                if lt_pot(pair2.bot_sig, reducer_sig)
                    reducer_sig = pair2.bot_sig
                    reducer_ind = pair2.bot_index
                end
            end
        end

        # add both as rows to matrix
        mult = divide(monomial(pair.top_sig),
                      monomial(basis.sigs[pair.top_index]))
        write_to_matrix_row!(matrix, basis, pair.top_index, symbol_ht,
                             ht, mult) 
        mult = divide(monomial(reducer_sig),
                      monomial(basis.sigs[reducer_ind]))
        lead_col_idx = write_to_matrix_row!(matrix, basis, reducer_ind,
                                            symbol_ht, ht, mult)

        # resize pivots array if needed
        if matrix.pivot_size < symbol_ht.load - 1
            resize!(matrix.pivots, 2 * (matrix.symbol_ht.load - 1))
            matrix.pivot_size = 2 * (matrix.symbol_ht.load - 1)
        end
        matrix.pivots[lead_col_idx] = matrix.nrows 
    end

    # remove selected pairs from pairset
    @inbounds for i in 1:pairset.load-npairs
        pairset.pairs[i] = pairset.pairs[i+npairs]
    end
    pairset.load -= npairs
end

function symbolic_pp!(basis::Basis{N},
                      matrix::MacaulayMatrix,
                      ht::MonomialHashtable,
                      symbol_ht::MonomialHashtable) where N

    pivots = matrix.pivots
    i = MonIdx(symbol_ht.offset)

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
            resize!(matrix.sig_order, matrix.size)
            resize!(matrix.coeffs, matrix.size)
        end

        exp = symbol_ht.exponents[i]
        divm = symbol_ht.hashdata[i].divmask
        mult = SVector{N, Exp}()
        
        j = 1
        @label target
        # find element in basis which divmask divides divmask of monomial
        @inbounds while j <= basis.basis_load && !divch(basis.lm_masks[j], divm)
            j += 1 
        end

        if j <= basis.basis_load
            @inbounds red_exp = leading_monomial(basis, ht, j)
            red_ind = j

            # actual divisibility check
            div_flag = true
            @inbounds for k in 1:N
                mult[k] = exp.exps[k] - red_exp.exps[k]
                if mult[k] < 0
                    div_flag = false
                    break
                end
            end
            if !div_flag
                j += 1
                @goto target
            end

            # set reducer
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
                mult2 = SVector{N, Exp}()
                @inbounds cand_exp = leading_monomial(basis, ht, j)
                @inbounds for k in 1:N
                    mult2[k] = exp.exps[k] - cand_exp.exps[k]
                    if mult2[k] < 0
                        div_flag = false
                        break
                    end
                end
                if !div_flag
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
                    comp_sigratio(basis, j, red_ind))
                    mult = mult2
                    red_ind = j
                    red_sig = cand_sig
                    @goto target2
                end
            end
            @inbounds lead_col_idx = write_to_matrix_row!(matrix, basis,
                                                          red_ind, symbol_ht,
                                                          ht, mult)
            pivots[lead_col_idx] = matrix.nrows
        end
    end
end

function write_to_matrix_row!(matrix::MacaulayMatrix,
                              basis::Basis,
                              basis_idx::Int,
                              symbol_ht::MonomialHashtable,
                              ht::MonomialHashtable,
                              mult::Monomial)

    hsh = Base.hash(mult)
    row_ind = matrix.nrows + 1
    poly = basis.monomials[basis_idx]
    row = similar(row)
    check_enlarge_hashtable!(symbol_ht, length(basis.monomials[basis_idx]))
    @inbounds matrix.rows[row_ind] =
        insert_multiplied_poly_in_hash_table!(row, hsh, mult, poly,
                                              ht, symbol_ht)
    matrix.coeffs[row_ind] = basis.coefficients[basis_idx]
    matrix.nrows += 1
    return first(matrix.rows[row_ind])
end