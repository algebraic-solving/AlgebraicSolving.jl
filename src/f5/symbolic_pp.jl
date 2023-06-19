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

function symbolic_pp!(basis::Basis{N},
                      matrix::MacaulayMatrix,
                      ht::MonomialHashtable,
                      symbol_ht::MonomialHashtable) where N

    pivots = matrix.pivots
    i = MonIdx(symbol_ht.offset)

    @inbounds while i <= symbol_ht.load
        if !iszero(matrix.pivots[i])
            i += one(MonIdx)
            continue
        end
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
        @inbounds while j <= basis.basis_load && !div(basis.lm_masks[j], divm)
            j += 1 
        end

        if j <= basis.basis_load
            @inbounds red_poly = basis.monomials[j]
            @inbounds red_sig = basis.sigs[j]
            @inbounds red_exp = ht.exponents[red_poly[1]]

            # actual divisibility check
            div_flag = true
            @inbounds for k in 1:N
                mult[k] = exp.exps[k] - red_exp.exps[k]
                if mult[k] < 0
                    div_flag = false
                    break
                end
            end
            j += 1
            if !div_flag
                @goto target
            end
            mul_red_sig = (index(red_sig), mul(mult, monomial(red_sig)))
            mul_sig_mask = divmask(monomial(mul_red_sig), ht.divmap,
                                   ht.ndivbits)

            # now we found a reducer
            j += 1
            @label target2
            @inbounds while j <= basis.basis_load && !div(basis.lm_masks[j], divm)
                j += 1 
            end

            mult2 = SVector{N, Exp}()
            @inbounds if j <= basis.basis_load
                # TODO first check if potential new reducer has smaller signature and divides lm
                # then do rewrite check
                cand_index = index(basis.sigmasks[j])
                cand_sig_mask = basis.sigmasks[j][2]
                if (cand_index == index(mul_red_sig)
                    && div(cand_sig_mask, mul_sig_mask))
                    cand_sig = basis.sigs[j]

                    m = monomial(mul_red_sig)
                    @inbounds for k in 1:N
                        mult2[k] = m.exps[k] - monomial(cand_sig).exps[k]
                        if mult2[k] < 0
                            j += 1
                            @goto target2
                        end
                    end
                    @inbounds red_poly = basis.monomials[j]
                    @inbounds red_exp = ht.exponents[red_poly[1]]
                    red_sig = cand_sig
                    mul_red_sig = (index(red_sig), mul(mult2, monomial(red_sig)))
                    mul_sig_mask = divmask(monomial(mul_red_sig), ht.divmap,
                                           ht.ndivbits)
                    j += 1
                    @goto target2
                end
            end
            # TODO: register red_poly, mul_red_sig etc
        end
    end
end
