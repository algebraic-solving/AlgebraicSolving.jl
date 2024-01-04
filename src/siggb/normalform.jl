function normalform(exps::Vector{Monomial{N}},
                    cfs::Vector{Coeff},
                    basis::Basis{N},
                    basis_ht::MonomialHashtable{N},
                    ind_order::IndOrder,
                    tags::Tags,
                    shift::Val{Shift},
                    vchar::Val{Char}) where {N, Char, Shift}

    symbol_ht = initialize_secondary_hash_table(basis_ht)
    matrix = initialize_matrix(Val(N))
    reinitialize_matrix!(matrix, 1)

    sig_ind = SigIndex(basis.input_load + 1)
    sig = (sig_ind, one_monomial(Monomial{N}))

    check_enlarge_hashtable!(symbol_ht, length(exps))
    @inbounds matrix.rows[1] = [insert_in_hash_table!(symbol_ht, e)
                                for e in exps]
    @inbounds matrix.coeffs[1] = cfs
    @inbounds matrix.sigs[1] = sig
    @inbounds matrix.parent_inds[1] = 0
    matrix.nrows = 1

    symbolic_pp!(basis, matrix, basis_ht, symbol_ht, ind_order, tags)
    @assert sig_ind == length(ind_order.ord) + 1
    push!(ind_order.ord, sig_ind)
    finalize_matrix!(matrix, symbol_ht, ind_order)
    pop!(ind_order.ord)

    echelonize!(matrix, tags, vchar, shift)

    # get result
    @inbounds res_exps = [symbol_ht.exponents[m_idx]
                for m_idx in matrix.rows[1]]
    @inbounds res_coeffs = matrix.coeffs[1]

    return res_exps, res_coeffs
end

# TODO: do we need to make sure that the product is normalized
function mult_pols(exps1::Vector{Monomial{N}},
                   exps2::Vector{Monomial{N}},
                   cfs1::Vector{Coeff},
                   cfs2::Vector{Coeff},
                   char::Val{Char}) where {N, Char}

    R, vrs = polynomial_ring(GF(Int(Char)), ["x$i" for i in 1:N],
                             ordering = :degrevlex)
    p1 = convert_to_pol(R, exps1, cfs1)
    p2 = convert_to_pol(R, exps2, cfs2)
    p = p1*p2

    lp = length(p)
    exps = exponent_vectors(p)
    cfs = coefficients(p)
    
    res_exps = Vector{Monomial{N}}(undef, lp)
    res_cfs = Vector{Coeff}(undef, lp)
    @inbounds for (i, (cf, evec)) in enumerate(zip(cfs, exps)) 
        m = monomial(SVector{N}((Exp).(evec)))
        cff = cf.data
        res_exps[i] = m
        res_cfs[i] = cff
    end

    return res_exps, res_cfs
end
