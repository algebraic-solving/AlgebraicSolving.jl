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

    echelonize!(matrix, tags, ind_order, vchar, shift, trace = false)

    # get result
    @inbounds res_exps = [symbol_ht.exponents[m_idx]
                for m_idx in matrix.rows[1]]
    @inbounds res_coeffs = matrix.coeffs[1]

    return res_exps, res_coeffs
end

# compute normal form with respect to basis
# *without* computing a GB for the corresponding ideal
function my_normal_form(mons::Vector{MonIdx},
                        coeffs::Vector{Coeff},
                        ht::MonomialHashtable{N},
                        gb::Vector{<:MPolyRingElem}) where N
    

    R = parent(first(gb))
    @inbounds mons = [ht.exponents[midx] for midx in mons]
    f = convert_to_pol(R, mons, coeffs)
    F = [f]
    return my_normal_form(F, gb)
end

function my_normal_form(mns::Vector{MonIdx},
                        cfs::Vector{Coeff},
                        bs_lens::Vector{Int32},
                        bs_cfs::Vector{Int32},
                        bs_exps::Vector{Int32},
                        basis_ht::MonomialHashtable{N},
                        char::Val{Char}) where {N, Char}

    nr_vars = N
    field_char = Int(Char)
    tbr_nr_gens = 1
    bs_nr_gens = length(bs_lens)
    is_gb = 1

    nr_thrds = 1
    info_level = 0

    mns_unhashed = [basis_ht.exponents[midx] for midx in mns]
    s = sortperm(mns_unhashed, lt = lt_drl, rev = true)
    mns_unhashed = mns_unhashed[s]
    cfs_s = copy(cfs)[s]
    tbr_lens, tbr_cfs, tbr_exps = _convert_to_msolve([mns_unhashed], [cfs_s])
                        
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
        nr_thrds, info_level)

    # convert to julia array, also give memory management to julia
    jl_ld   = nf_ld[]
    jl_len  = Base.unsafe_wrap(Array, nf_len[], jl_ld)
    jl_exp  = Base.unsafe_wrap(Array, nf_exp[], nr_terms*nr_vars)
    ptr     = reinterpret(Ptr{Int32}, nf_cf[])
    jl_cf   = Base.unsafe_wrap(Array, ptr, nr_terms)

    coeffs_res, mns_res = convert_ms_to_ht(jl_ld, jl_len, jl_cf,
                                           jl_exp,
                                           basis_ht)

    ccall((:free_f4_julia_result_data, libneogb), Nothing ,
          (Ptr{Nothing}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}},
           Ptr{Ptr{Cvoid}}, Int, Int),
          cglobal(:jl_free), nf_len, nf_exp, nf_cf, jl_ld, field_char)

    return first(coeffs_res), first(mns_res)
end

function my_normal_form(F::Vector{T},
                        gb::Vector{T};
                        nr_thrds::Int=1,
                        info_level::Int=0) where T <: MPolyRingElem

    if (length(F) == 0 || length(gb) == 0)
        error("Input data not valid.")
    end

    R = first(F).parent
    @assert Nemo.ordering(R) == :degrevlex
    nr_vars     = nvars(R)
    field_char  = Int(characteristic(R))

    if !(is_probable_prime(field_char))
        error("At the moment we only supports finite fields.")
    end

    G = gb

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
        nr_thrds, info_level)

    # convert to julia array, also give memory management to julia
    jl_ld   = nf_ld[]
    jl_len  = Base.unsafe_wrap(Array, nf_len[], jl_ld)
    jl_exp  = Base.unsafe_wrap(Array, nf_exp[], nr_terms*nr_vars)
    ptr     = reinterpret(Ptr{Int32}, nf_cf[])
    jl_cf   = Base.unsafe_wrap(Array, ptr, nr_terms)

    basis = _convert_finite_field_array_to_abstract_algebra(
                jl_ld, jl_len, jl_cf, jl_exp, R, 0)

    ccall((:free_f4_julia_result_data, libneogb), Nothing ,
          (Ptr{Nothing}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}},
           Ptr{Ptr{Cvoid}}, Int, Int),
          cglobal(:jl_free), nf_len, nf_exp, nf_cf, jl_ld, field_char)

    return basis
end

# compute normal forms with msolve
@inline function _convert_to_msolve(exps::Vector{Vector{T}},
                                    cfs::Vector{Vector{Coeff}}) where {T <: Monomial}

    len = length(exps)
    @inbounds ms_cfs = (Int32).(vcat(cfs...))
    @inbounds ms_exps = vcat([vcat([convert(Vector{Int32}, ex.exps) for ex in exps[i]]...)
                              for i in 1:length(exps)]...)
    return (Int32).([length(ex) for ex in exps]), ms_cfs, ms_exps
end

function convert_ms_to_ht(bld::Int32,
                          blen::Vector{Int32},
                          bcf::Vector{Int32},
                          bexp::Vector{Int32},
                          ht::MonomialHashtable{N}) where N

    nr_gens = bld

    coeffs_res = Vector{Vector{Coeff}}(undef, nr_gens)
    mns_res = Vector{Vector{MonIdx}}(undef, nr_gens)

    len = 0
    for i in 1:nr_gens
        cfs_i = Vector{Coeff}(undef, blen[i])
        mns_i = Vector{MonIdx}(undef, blen[i])
        coeffs_res[i] = cfs_i
        mns_res[i] = mns_i
        if bcf[len+1] == 0
            resize!(cfs_i, 0)
            resize!(mns_i, 0)
            continue
        else
            for j in 1:blen[i]
                cfs_i[j] = Coeff(bcf[len+j]) 
                @assert !iszero(cfs_i[j])
                exp_v = SVector{N, Exp}(bexp[(len+j-1)*N+1:(len+j)*N])
                mns_i[j] = insert_in_hash_table!(ht, monomial(exp_v))
            end
        end
        len += blen[i]
    end

    return coeffs_res, mns_res
end

function _convert_basis_to_msolve(basis::Basis,
                                  basis_ht::MonomialHashtable)

    l = basis.basis_load - basis.basis_offset + 1
    lens = Vector{Int32}(undef, l)
    ms_cfs = Int32[]
    ms_exps = Int32[]

    @inbounds for i in basis.basis_offset:basis.basis_load
        exps = [basis_ht.exponents[eidx] for eidx in basis.monomials[i]]
        cfs = basis.coefficients[i]
        plen, pms_cfs, pms_exps = _convert_to_msolve(exps, cfs)
        append!(lens, plen)
        append!(ms_cfs, pms_cfs)
        append!(ms_exps, pms_exps)
    end

    return lens, ms_cfs, ms_exps
end
