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

    symbolic_pp!(basis, matrix, basis_ht, symbol_ht, ind_order, tags,
                 forbidden_tag = :col)

    push!(ind_order.ord, ind_order.max_ind + one(SigIndex))
    finalize_matrix!(matrix, symbol_ht, ind_order)
    pop!(ind_order.ord)

    echelonize!(matrix, tags, vchar, shift, trace = false)

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

# compute normal form with respect to basis
# *without* computing a GB for the corresponding ideal
function my_iszero_normal_form(mons::Vector{Monomial{N}},
                               coeffs::Vector{Coeff},
                               basis::Basis,
                               basis_ht::MonomialHashtable,
                               char::Val{Char}) where {N, Char}
    

    R, vrs = polynomial_ring(GF(Int(Char)), ["x$i" for i in 1:N],
                             ordering = :degrevlex)
    f = convert_to_pol(R, mons, coeffs)
    F = [f]
    nr_vars     = nvars(R)
    field_char  = Int(characteristic(R))

    #= first get a degree reverse lexicographical Gröbner basis for I =#
    G0 = [convert_to_pol(R, [basis_ht.exponents[midx] for midx in basis.monomials[i]],
                         basis.coefficients[i]) for i in basis.basis_offset:basis.basis_load]
    G = groebner_basis(Ideal(G0), complete_reduction = true)
    # G = G0

    tbr_nr_gens = length(F)
    bs_nr_gens  = length(G)
    is_gb       = 1

    # convert ideal to flattened arrays of ints
    tbr_lens, tbr_cfs, tbr_exps = _convert_to_msolve(F)
    bs_lens, bs_cfs, bs_exps    = _convert_to_msolve(G)

    nf_ld  = Ref(Cint(0))
    nf_len = Ref(Ptr{Cint}(0))
    nf_exp = Ref(Ptr{Cint}(0))
    nf_cf  = Ref(Ptr{Cvoid}(0))

    nr_terms  = ccall((:export_nf, libneogb), Int,
        (Ptr{Nothing}, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cvoid}},
        Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cvoid},
        Cint, Cint, Cint, Cint, Cint, Cint, Cint),
        cglobal(:jl_malloc), nf_ld, nf_len, nf_exp, nf_cf, tbr_nr_gens, tbr_lens, tbr_exps,
        tbr_cfs, bs_nr_gens, bs_lens, bs_exps, bs_cfs, field_char, 0, 0, nr_vars, is_gb,
        1, 0)

    # convert to julia array, also give memory management to julia
    jl_ld   = nf_ld[]
    jl_len  = Base.unsafe_wrap(Array, nf_len[], jl_ld)
    jl_exp  = Base.unsafe_wrap(Array, nf_exp[], nr_terms*nr_vars)
    ptr     = reinterpret(Ptr{Int32}, nf_cf[])
    jl_cf   = Base.unsafe_wrap(Array, ptr, nr_terms)

    res = _convert_finite_field_array_to_abstract_algebra(
        jl_ld, jl_len, jl_cf, jl_exp, R, 0)

    ccall((:free_f4_julia_result_data, libneogb), Nothing ,
          (Ptr{Nothing}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}},
           Ptr{Ptr{Cvoid}}, Int, Int),
          cglobal(:jl_free), nf_len, nf_exp, nf_cf, jl_ld, field_char)

    return iszero(first(res))
end
