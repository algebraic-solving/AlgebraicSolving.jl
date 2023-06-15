# TODO: enforce some type restrictions
# TODO: make sure pivots array has enough space
function select_normal!(pairset::Pairset,
                        basis::Basis,
                        matrix::MacaulayMatrix,
                        ht::MonomialHashtable,
                        symbol_ht::MonomialHashtable) 

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

    reinitialize_matrix!(matrix, npairs)
    skip = falses(npairs)

    @inbounds for i in 1:npairs
        skip[i] && continue
        pair = pairset.pairs[i]
        curr_top_sig = pair.top_sig
        reducer_sig = pair.bot_sig
        reducer_ind = pair.bot_index

        @inbounds for j in (i+1):npairs
            pair2 = pairset.pairs[j]
            if pair2.top_sig = curr_top_sig
                skip[j] = true
                if lt_pot(pair2.bot_sig, reducer_sig)
                    reducer_sig = pair2.bot_sig
                    reducer_ind = pair2.bot_index
                end
            end
        end
        write_to_matrix_row!(matrix, basis, pair.top_index, symbol_ht,
                             ht, curr_top_sig[2]) 
        lead_col_idx = write_to_matrix_row!(matrix, basis, reducer_ind,
                                            symbol_ht,
                                            ht, reducer_sig[2])
        pivots[lead_col_idx] = matrix.nrows 
    end

    # remove selected parirs from pairset
    @inbounds for i in 1:pairset.load-npairs
        pairset.pairs[i] = pairset.[i+npairs]
    end
    pairset.load -= npairs
end

# TODO: make sure to have space before calling this
function write_to_matrix_row!(matrix::MacaulayMatrix,
                              basis::Basis,
                              basis_idx::Int,
                              symbol_ht::MonomialHashtable,
                              ht::MonomialHashtable,
                              top_sig_mon::Monomial)

    mult = divide(top_sig_mon, basis.sigs[basis_idx][2])
    hsh = Base.hash(mult)
    row_ind = matrix.nrows + 1
    @inbounds matrix.rows[row_ind] =
        multiplied_poly_to_matrix_row!(symbol_ht, ht,
                                       hsh, mult, basis.monomials[basis_idx])
    matrix.coeffs[row_ind] = basis.coefficients[basis_idx]
    matrix.nrows += 1
    return first(matrix.rows[row_ind])
end
