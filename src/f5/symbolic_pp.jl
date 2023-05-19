function select_normal!(pairset::Pairset,
                        basis::POTBasis,
                        matrix::MacaulayMatrix,
                        ht::MonomialHashtable,
                        symbol_ht::MonomialHashtable) 

    # number of selected pairs
    npairs = 0
    sig_deg = zero(Exp) 
    for i in 1:pairset.load
        if iszero(sig_deg)
            sig_deg = pairset.pairs[i].top_sig.deg
            npairs += 1
            iszero(sig_deg) && break
        end
        if pairset.pairs[i].top_sig.deg == deg
            npairs += 1
        else
            break
        end
    end

    reinitialize_matrix!(matrix, npairs)

    # TODO: EVERYTHING BELOW HERE STILL FROM GROEBNER.JL
    # polynomials from pairs in order (p11, p12)(p21, p21)
    # (future rows of the matrix)
    gens = Vector{Int}(undef, 2 * npairs)

    # monomial buffer
    etmp = ht.exponents[1]
    i = 1
    @inbounds while i <= npairs
        matrix.ncols += 1
        load = 1
        lcm = ps[i].lcm
        j = i

        # we collect all generators with same lcm into gens
        while j <= npairs && ps[j].lcm == lcm
            gens[load] = ps[j].poly1
            load += 1
            gens[load] = ps[j].poly2
            load += 1
            j += 1
        end
        load -= 1

        # sort by the index in the basis (by=identity)
        sort_generators_by_position!(gens, load)

        # now we collect reducers and to-be-reduced polynomials

        # first generator index in groebner basis
        prev = gens[1]
        # first generator in hash table
        poly = basis.monoms[prev]
        # first generator lead monomial index in hash data
        vidx = poly[1]

        # first generator exponent
        eidx = ht.exponents[vidx]
        # exponent of lcm corresponding to first generator
        elcm = ht.exponents[lcm]
        etmp = monom_division!(etmp, elcm, eidx)
        # now etmp contents complement to eidx in elcm

        # hash of complement
        htmp = ht.hashdata[lcm].hash - ht.hashdata[vidx].hash

        # add row as a reducer
        matrix.nup += 1
        uprows[matrix.nup] = multiplied_poly_to_matrix_row!(symbol_ht, ht, htmp, etmp, poly)
        # map upper row to index in basis
        matrix.up2coef[matrix.nup] = prev

        # mark lcm column as reducer in symbolic hashtable
        symbol_ht.hashdata[uprows[matrix.nup][1]].idx = 2
        # increase number of rows set
        matrix.nrows += 1

        # over all polys with same lcm,
        # add them to the lower part of matrix
        @inbounds for k in 1:load
            # duplicate generator,
            # we can do so as long as generators are sorted
            if gens[k] == prev
                continue
            end

            # if the table was reallocated
            elcm = ht.exponents[lcm]

            # index in gb
            prev = gens[k]
            # poly of indices of monoms in hash table
            poly = basis.monoms[prev]
            vidx = poly[1]
            # leading monom idx
            eidx = ht.exponents[vidx]

            etmp = monom_division!(etmp, elcm, eidx)

            htmp = ht.hashdata[lcm].hash - ht.hashdata[vidx].hash

            # add row to be reduced
            matrix.nlow += 1
            lowrows[matrix.nlow] = multiplied_poly_to_matrix_row!(symbol_ht, ht, htmp, etmp, poly)
            # map lower row to index in basis
            matrix.low2coef[matrix.nlow] = prev

            symbol_ht.hashdata[lowrows[matrix.nlow][1]].idx = 2

            matrix.nrows += 1
        end

        i = j
    end

    resize!(matrix.lowrows, matrix.nrows - matrix.ncols)

    # remove selected parirs from pairset
    @inbounds for i in 1:pairset.load-npairs
        ps[i] = ps[i+npairs]
    end
    pairset.load -= npairs
end
