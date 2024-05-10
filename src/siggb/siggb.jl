# TODO:
# insert elems in nondeg one by one
# divmask at random two times (for example)?


# sizes for initialization
const init_ht_size = 2^17
const init_basis_size = 10000
const init_syz_size = 1000
const init_pair_size = 10000
# default sorting alg
const def_sort_alg = Base.Sort.DEFAULT_UNSTABLE
include("typedefs.jl")
include("monomials.jl")
include("hashtable.jl")
include("rewriting.jl")
include("update.jl")
include("symbolic_pp.jl")
include("linear_algebra.jl")
include("tracer.jl")
include("module.jl")
include("normalform.jl")
include("affine_cells.jl")
include("interfaces.jl")


#---------------- user functions --------------------#

@doc Markdown.doc"""
function sig_groebner_basis(sys::Vector{T}; info_level::Int=0, degbound::Int=0, mod_ord::Symbol=:DPOT) where {T <: MPolyRingElem}

Compute a Signature Gröbner basis of the sequence `sys` w.r.t. to the
degree reverse lexicographical monomial ordering and the degree
position-over-term ordering induced by `sys`. The output is a vector
of `Tuple{Tuple{Int64, T}, T}` where the first element indicates the
signature and the second the underlying polynomial.

**Note**: At the moment only ground fields of characteristic `p`, `p` prime, `p < 2^{31}` are supported.
**Note**: The input generators must be homogeneous.
**Note**: The algorithms behaviour may depend heavily on how the elements in `sys` are sorted.

# Arguments
- `sys::Vector{T} where T <: MpolyElem`: input generators.
- `info_level::Int=0`: info level printout: off (`0`, default), computational details (`1`)
- `degbound::Int=0`: Compute a full Gröbner basis if `0` otherwise only up to degree `degbound`.
- `mod_ord::Symbol=:DPOT`: The module monomial order underlying the computation, either `:DPOT` (default) or `:POT`.

# Example
```jldoctest
julia> using AlgebraicSolving

julia> R, vars = polynomial_ring(GF(17), ["x$i" for i in 1:4])
(Multivariate polynomial ring in 4 variables over GF(17), fpMPolyRingElem[x1, x2, x3, x4])

julia> F = AlgebraicSolving.cyclic(R)
fpMPolyRingElem[x1 + x2 + x3 + x4, x1*x2 + x1*x4 + x2*x3 + x3*x4, x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4, x1*x2*x3*x4 + 16]

julia> Fhom = AlgebraicSolving._homogenize(F.gens)
4-element Vector{fpMPolyRingElem}:
 x1 + x2 + x3 + x4
 x1*x2 + x2*x3 + x1*x4 + x3*x4
 x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4
 x1*x2*x3*x4 + 16*x5^4

julia> sig_groebner_basis(Fhom)
7-element Vector{Tuple{Tuple{Int64, fpMPolyRingElem}, fpMPolyRingElem}}:
 ((1, 1), x1 + x2 + x3 + x4)
 ((2, 1), x2^2 + 2*x2*x4 + x4^2)
 ((3, 1), x2*x3^2 + x3^2*x4 + 16*x2*x4^2 + 16*x4^3)
 ((4, 1), x2*x3*x4^2 + x3^2*x4^2 + 16*x2*x4^3 + x3*x4^3 + 16*x4^4 + 16*x5^4)
 ((4, x3), x3^3*x4^2 + x3^2*x4^3 + 16*x3*x5^4 + 16*x4*x5^4)
 ((4, x2), x2*x4^4 + x4^5 + 16*x2*x5^4 + 16*x4*x5^4)
 ((4, x2*x3), x3^2*x4^4 + x2*x3*x5^4 + 16*x2*x4*x5^4 + x3*x4*x5^4 + 15*x4^2*x5^4)
```
"""
function sig_groebner_basis(sys::Vector{T}; info_level::Int=0, degbound::Int=0, mod_ord::Symbol=:DPOT) where {T <: MPolyRingElem}

    # data structure setup/conversion
    sys_mons, sys_coeffs, basis_ht, char, shift = input_setup(sys, mod_ord)

    # fill basis, pairset, tags
    basis, pairset, tags, ind_order, tr = fill_structs!(sys_mons, sys_coeffs, basis_ht)

    sysl = length(sys)
    # compute divmasks
    fill_divmask!(basis_ht)
    @inbounds for i in 1:sysl
        basis.lm_masks[i] = basis_ht.hashdata[basis.monomials[i][1]].divmask
    end

    logger = ConsoleLogger(stdout, info_level == 0 ? Warn : Info)
    with_logger(logger) do
        _, arit_ops = siggb!(basis, pairset, basis_ht, char, shift,
                             tags, ind_order, tr, degbound,
                             mod_ord)
        @info "$(arit_ops) total submul's"
    end

    # output
    R = parent(first(sys))
    eltp = typeof(first(sys))
    outp = Tuple{Tuple{Int, eltp}, eltp}[]
    @inbounds for i in basis.basis_offset:basis.basis_load
        pol = convert_to_pol(R,
                             [basis_ht.exponents[m] for m in basis.monomials[i]],
                             basis.coefficients[i])
        s = basis.sigs[i]
        ctx = MPolyBuildCtx(R)
        push_term!(ctx, base_ring(R)(1), Vector{Int}(monomial(s).exps))
        sig = (Int(index(s)), finish(ctx))

        push!(outp, (sig, pol))
    end

    return outp
end

function sig_sat(sys::Vector{T}, H::Vector{T}; info_level::Int=0) where {T <: MPolyRingElem}

    full_sys = vcat(sys, H)
    
    # data structure setup/conversion
    sys_mons, sys_coeffs, basis_ht, char, shift = input_setup(full_sys)

    # fill basis, pairset, tags
    basis, pairset, tags, ind_order, tr = fill_structs!(sys_mons, sys_coeffs, basis_ht)

    # set tag for sats
    sysl = length(full_sys)
    for idx in length(sys)+1:sysl
        tags[idx] = :sat
    end
    make_sat_incompat!(tags, ind_order)

    # compute divmasks
    fill_divmask!(basis_ht)
    @inbounds for i in 1:sysl
        basis.lm_masks[i] = basis_ht.hashdata[basis.monomials[i][1]].divmask
    end

    logger = ConsoleLogger(stdout, info_level == 0 ? Warn : Info)
    with_logger(logger) do
        added_unit, _ = siggb!(basis, pairset, basis_ht, char, shift,
                               tags, ind_order, tr, trace=true)
        # output
        R = parent(first(sys))
        eltp = typeof(first(sys))
        outp = eltp[]
        if added_unit
            push!(outp, one(R))
        else
            @inbounds for i in basis.basis_offset:basis.basis_load
                gettag(tags, index(basis.sigs[i])) == :sat && continue
                pol = convert_to_pol(R,
                                     [basis_ht.exponents[m] for m in basis.monomials[i]],
                                     basis.coefficients[i])
                push!(outp, pol)
            end
        end
        return outp
    end
end

function nondeg_locus(sys::Vector{T}; info_level::Int=0) where {T <: MPolyRingElem}

    # data structure setup/conversion
    sys_mons, sys_coeffs, basis_ht, char, shift = input_setup(sys, :POT)

    # fill basis, pairset, tags
    basis, _, tags, ind_order, tr = fill_structs!(sys_mons, sys_coeffs, basis_ht, def_tg=:ndeg)
    R = parent(first(sys))
    nvrs = ngens(R)
    r_char = characteristic(R)
    pairset = init_pairset(Val(nvrs))

    # compute divmasks
    fill_divmask!(basis_ht)
    sysl = length(sys)
    @inbounds for i in 1:sysl
        basis.lm_masks[i] = basis_ht.hashdata[basis.monomials[i][1]].divmask
    end

    logger = ConsoleLogger(stdout, info_level == 0 ? Warn : Info)
    with_logger(logger) do
        arit_ops = 0
        ind_conn = IndConn()
        for i in 1:sysl
            @info "sequence index $(i)"
            add_unit_pair!(basis, pairset, i, basis.degs[i])
            ndeg_ins_index = i < sysl ? SigIndex(i+1) : zero(SigIndex)
            added_unit, new_arit_ops = siggb!(basis, pairset, basis_ht, char, shift,
                                              tags, ind_order, tr, 0, :POT, true,
                                              ndeg_ins_index, ind_conn)
            arit_ops += new_arit_ops

            # extract suitable random linear combinations of inserted elements
            sig_inds = (SigIndex).(collect(eachindex(ind_order.ord)))
            for idx in sig_inds
                gettag(tags, idx) != :fndegins && continue
                # check if element was inserted in previous round
                gr_inds = filter(idx2 -> gettag(tags, idx2) == :ndeg && !cmp_ind(idx2, idx, ind_order),
                                 sig_inds)
                f_ord_ind = ind_order.ord[i]
                min_gr_ord_ind, _ = findmin(idx2 -> ind_order.ord[idx2], gr_inds)
                f_ord_ind != min_gr_ord_ind && continue

                # construct nonzero condition to append
                nz_cfs, nz_mons = Coeff[], MonIdx[]
                for (j, syz_msk) in enumerate(basis.syz_masks[1:basis.syz_load])
                    if index(syz_msk) == idx
                        to_add_cfs, to_add_mns = construct_module((idx, basis.syz_sigs[j]),
                                                                  basis,
                                                                  basis_ht,
                                                                  tr.syz_ind_to_mat[j],
                                                                  tr,
                                                                  char,
                                                                  ind_order, idx)
                        mul_by_coeff!(to_add_cfs, rand(one(Coeff):Coeff(r_char-1)),
                                      char)
                        nz_cfs, nz_mons = add_pols(nz_cfs, nz_mons,
                                                   to_add_cfs, to_add_mns, char)
                    end
                end
                is_one((nz_cfs, nz_mons), basis_ht) && continue
                sort_poly!((nz_cfs, nz_mons), by = midx -> basis_ht.exponents[midx],
                           lt = lt_drl, rev = true)
                normalize_cfs!(nz_cfs, char)
                add_new_sequence_element!(basis, basis_ht, tr, nz_cfs, nz_mons,
                                          ind_order, ind_order.max_ind+one(SigIndex),
                                          pairset, tags,
                                          new_tg = :sat)
                make_sat_incompat!(tags, ind_order)
            end
            @info "------------------------------------------"
            @info "saturation step"
            _, new_arit_ops = siggb!(basis, pairset, basis_ht, char, shift,
                                     tags, ind_order, tr, 0, :POT, true,
                                     ndeg_ins_index, ind_conn)
            arit_ops += new_arit_ops
            @info "------------------------------------------"
            @info "$(arit_ops) total submul's"
        end

        # output
        R = parent(first(sys))
        eltp = typeof(first(sys))
        outp = eltp[]
        @inbounds for i in basis.basis_offset:basis.basis_load
            gettag(tags, index(basis.sigs[i])) == :sat && continue
            pol = convert_to_pol(R,
                                 [basis_ht.exponents[m] for m in basis.monomials[i]],
                                 basis.coefficients[i])
            push!(outp, pol)
        end
        return outp
    end
end

function sig_decomp(sys::Vector{T}; info_level::Int=0) where {T <: MPolyRingElem}

    # data structure setup/conversion
    sys_mons, sys_coeffs, basis_ht, char, shift = input_setup(sys)
    
    # fill basis, pairset, tags
    basis, pairset, tags, ind_order, tr = fill_structs!(sys_mons, sys_coeffs,
                                                        basis_ht, def_tg=:split)

    sysl = length(sys)
    # compute divmasks
    fill_divmask!(basis_ht)
    @inbounds for i in 1:sysl
        basis.lm_masks[i] = basis_ht.hashdata[basis.monomials[i][1]].divmask
    end

    logger = ConsoleLogger(stdout, info_level == 0 ? Warn : Info)
    result = with_logger(logger) do
        R = parent(first(sys))
        timer = new_timer()
        lc_sets = sig_decomp!(basis, pairset, basis_ht, char, shift,
                              tags, ind_order, tr, R, timer)
        @info timer
        return lc_sets
    end
end

#---------------- function for sig_groebner_basis --------------------#

function siggb!(basis::Basis{N},
                pairset::Pairset,
                basis_ht::MonomialHashtable,
                char::Val{Char},
                shift::Val{Shift},
                tags::Tags,
                ind_order::IndOrder,
                tr::Tracer,
                degbound::Int=0,
                mod_ord::Symbol=:DPOT,
                ndeg_ins_index::SigIndex=zero(SigIndex),
                conn_indices::IndConn=IndConn()) where {N, Char, Shift}

    # fake syz queue
    syz_queue = SyzInfo[]
    arit_ops = 0

    sort_pairset!(pairset, 1, pairset.load-1, mod_ord, ind_order)

    while !iszero(pairset.load)
        if !iszero(degbound) && first(pairset.elems).deg > degbound
            break
        end
	matrix = initialize_matrix(Val(N))
        symbol_ht = initialize_secondary_hash_table(basis_ht)

        _, compat_ind, sigind = select_normal!(pairset, basis, matrix,
                                               basis_ht, symbol_ht, ind_order, tags,
                                               mod_ord)
        symbolic_pp!(basis, matrix, basis_ht, symbol_ht,
                     ind_order, tags, sigind, compat_ind,
                     mod_ord)
        finalize_matrix!(matrix, symbol_ht, ind_order)
        iszero(matrix.nrows) && continue
        arit_ops_new = echelonize!(matrix, tags, ind_order, char,
                                   shift, trace)
        arit_ops += arit_ops_new

        added_unit = update_siggb!(basis, matrix, pairset, symbol_ht,
                                   basis_ht, ind_order, tags,
                                   tr, char, syz_queue, conn_indices,
                                   trace,
                                   ndeg_ins_index)
        if added_unit
            return true, arit_ops
        end
        sort_pairset!(pairset, 1, pairset.load-1, mod_ord, ind_order)
    end

    return false, arit_ops
end

#---------------- functions for splitting --------------------#

function sig_decomp!(basis::Basis{N},
                     pairset::Pairset,
                     basis_ht::MonomialHashtable,
                     char::Val{Char},
                     shift::Val{Shift},
                     tags::Tags,
                     ind_order::IndOrder,
                     tr::SigTracer,
                     R::MPolyRing,
                     timer::Timings) where {N, Char, Shift}

    # compute ideal
    eqns = [convert_to_pol(R, [basis_ht.exponents[mdx] for mdx in basis.monomials[i]],
                           basis.coefficients[i])
            for i in 1:basis.input_load]
    X = WLocClosedSet{eltype(eqns)}(eqns)
    
    queue = [(basis, pairset, tags, ind_order, X, SyzInfo[], tr)]
    result = WLocClosedSet[]

    while !isempty(queue)
        bs, ps, tgs, ind_ord, lc_set, syz_queue, tr = popfirst!(queue)
        neqns = num_eqns(lc_set)
        filter!(gb -> !(one(R) in gb), lc_set.gbs)
        @info "starting component, $(length(queue)) remaining, $(neqns) equations"
        if any(isone, lc_set.hull_eqns)
            @info "is empty set"
            @info "------------------------------------------"
            continue
        elseif add_to_output!(result, lc_set)
            @info "------------------------------------------"
            continue
        end
        found_zd, isempt, zd_coeffs,
        zd_mons, zd_ind = siggb_for_split!(bs, ps,
                                           tgs, ind_ord,
                                           basis_ht, tr,
                                           syz_queue,
                                           char, shift, lc_set,
                                           timer)
        if found_zd
            @info "splitting component"
            tim = @elapsed lc_set_hull, bs2, ps2, tgs2,
                           ind_ord2, lc_set_nz, tr2 = split!(bs, basis_ht,
                                                             zd_mons, zd_coeffs,
                                                             tr, ps,
                                                             zd_ind, tgs,
                                                             ind_ord,
                                                             lc_set)
            timer.comp_lc_time += tim
            pushfirst!(queue, (bs, ps, tgs, ind_ord, lc_set_hull, syz_queue, tr))
            pushfirst!(queue, (bs2, ps2, tgs2, ind_ord2, lc_set_nz, SyzInfo[], tr2))
        else
            @info "finished component"
            push!(result, lc_set)
        end
        @info "------------------------------------------"
    end
    return result
end

function siggb_for_split!(basis::Basis{N},
                          pairset::Pairset,
                          tags::Tags,
                          ind_order::IndOrder,
                          basis_ht::MonomialHashtable,
                          tr::SigTracer,
                          syz_queue::Vector{SyzInfo},
                          char::Val{Char},
                          shift::Val{Shift},
                          lc_set::AffineCell,
                          timer::Timings) where {N, Char, Shift}

    splitting_inds = [index(basis.sigs[i]) for i in 1:basis.input_load]
    filter!(ind -> gettag(tags, ind) == :split, splitting_inds)
    sort!(splitting_inds, by = ind -> ind_order.ord[ind])

    sort_pairset!(pairset, 1, pairset.load-1, ind_order)

    lt_squeue = (syz1, syz2) -> begin
        idx1 = index(syz1)
        idx2 = index(syz2)
        mon1 = monomial(syz1)
        mon2 = monomial(syz2)
        if mon1.deg == mon2.deg
            ind_order.ord[idx1] < ind_order.ord[idx2]
        else
            mon1.deg < mon2.deg
        end
    end

    while !iszero(pairset.load)
        # find minimum pair index
        min_pair_idx = minimum(pair -> ind_order.ord[index(pair.top_sig)],
                               pairset.elems[1:pairset.load])

	matrix = initialize_matrix(Val(N))
        symbol_ht = initialize_secondary_hash_table(basis_ht)

        tim = @elapsed deg, _ = select_normal!(pairset, basis, matrix,
                                               basis_ht, symbol_ht, ind_order, tags)
        timer.select_time += tim
        tim = @elapsed symbolic_pp!(basis, matrix, basis_ht, symbol_ht,
                                    ind_order, tags)
        timer.sym_pp_time += tim

        finalize_matrix!(matrix, symbol_ht, ind_order)
        iszero(matrix.nrows) && continue
        tim = @elapsed echelonize!(matrix, tags, ind_order, char, shift)
        timer.lin_alg_time += tim

        time = @elapsed update_siggb!(basis, matrix, pairset,
                                      symbol_ht, basis_ht,
                                      ind_order, tags,
                                      tr, char, syz_queue)
        timer.update_time += tim

        # find minimum syzygy index
        sort!(syz_queue, by = sz -> (basis.syz_sigs[sz[1]].deg, ind_order.ord[index(basis.syz_masks[sz[1]])]))
        if !isempty(syz_queue)
            min_syz_idx = minimum(sz -> ind_order.ord[index(basis.syz_masks[sz[1]])], syz_queue)
        else
            min_syz_idx = min_pair_idx
        end
        regular_up_to = min(min_pair_idx, min_syz_idx) - 1

        # check to see if we can split with one of the syzygies
        if !isempty(syz_queue)
            fs_syz_idx = first(syz_queue)[1]
            min_syz_deg = basis.syz_sigs[fs_syz_idx].deg
            poss_syz_new_deg = deg - basis.degs[basis.input_load]
            if min_syz_deg <= poss_syz_new_deg
                does_split, cofac_coeffs,
                cofac_mons, cofac_ind = process_syz_for_split!(syz_queue, basis_ht,
                                                               basis, tr, ind_order, char, lc_set,
                                                               tags, splitting_inds, regular_up_to,
                                                               timer)
                if does_split
                    return true, false, cofac_coeffs,
                    cofac_mons, cofac_ind
                end
            end
        end

        sort_pairset!(pairset, 1, pairset.load-1, ind_order)
    end
    if !isempty(syz_queue)
        sort!(syz_queue, by = sz -> basis.syz_sigs[sz[1]].deg)
        regular_up_to = minimum(sz -> ind_order.ord[index(basis.syz_masks[sz[1]])], syz_queue) - 1
        does_split, cofac_coeffs, cofac_mons,
        cofac_ind = process_syz_for_split!(syz_queue, basis_ht,
                                           basis, tr, ind_order, char, lc_set,
                                           tags, splitting_inds, regular_up_to,
                                           timer)
        if does_split
            return true, false, cofac_coeffs,
            cofac_mons, cofac_ind
        end
    end

    return false, false, Coeff[], Monomial{N}[], zero(SigIndex)
end

function split!(basis::Basis{N},
                basis_ht::MonomialHashtable{N},
                cofac_mons_hsh::Vector{MonIdx},
                cofac_coeffs::Vector{Coeff},
                tr::SigTracer,
                pairset::Pairset{N},
                zd_ind::SigIndex,
                tags::Tags,
                ind_order::IndOrder,
                lc_set::AffineCell) where N


    @inbounds begin
        # new polynomial
        cofac_mons = [basis_ht.exponents[midx] for midx in cofac_mons_hsh]
        h = convert_to_pol(ring(lc_set), cofac_mons, cofac_coeffs)
        @info "splitting index $(zd_ind), degree $(total_degree(h))"

        # component with nonzero condition
        sorted_inds = collect(1:basis.input_load)
        to_del = findall(idx -> basis.is_red[idx] || gettag(tags, idx) == :hull, sorted_inds)
        push!(to_del, zd_ind)
        unique!(sort!(to_del))
        deleteat!(sorted_inds, to_del)
        sort!(sorted_inds, by = ind -> ind_order.ord[ind])

        sys2_mons = copy(basis.monomials[sorted_inds])
        sys2_coeffs = copy(basis.coefficients[sorted_inds])
        basis2, ps2, tags2, ind_ord2, tr2 = fill_structs!(sys2_mons, sys2_coeffs,
                                                          basis_ht, def_tg=:split)

        ind_ord2 = new_ind_order(basis2)
        
        # hull component
        zd_deg = basis_ht.exponents[first(cofac_mons_hsh)].deg

        ge_deg_inds = filter(i -> basis.degs[i] > zd_deg, 1:basis.input_load)
        if isempty(ge_deg_inds)
            ord_ind = ind_order.max_ind + one(SigIndex)
        else
            ord_ind, _ = findmin(i -> ind_order.ord[i], ge_deg_inds)
        end

        # insert zd in system
        s_ind = add_new_sequence_element!(basis, basis_ht, tr,
                                          cofac_coeffs, cofac_mons_hsh,
                                          ind_order, ord_ind, pairset,
                                          tags, new_tg = :split)

        # new components
        lc_set_hull, lc_set_nz = split(lc_set, h)
        push!(lc_set_hull.seq, h)
        push!(lc_set_hull.tags, :hull)
        push!(lc_set_hull.hull_eqns, h)

        # TODO: do we need to delete from hull eqns?
        deleteat!(lc_set_nz.seq, to_del)
        deleteat!(lc_set_nz.tags, to_del)
        new_codim_ub = min(codim_upper_bound(lc_set), num_eqns(lc_set) - 1)
        add_hyperplanes!(lc_set_nz, new_codim_ub)
    end

    return lc_set_hull, basis2, ps2, tags2, ind_ord2, lc_set_nz, tr2
end    

function process_syz_for_split!(syz_queue::Vector{SyzInfo},
                                basis_ht::MonomialHashtable,
                                basis::Basis{N},
                                tr::SigTracer,
                                ind_order::IndOrder,
                                char::Val{Char},
                                lc_set::AffineCell,
                                tags::Tags,
                                splitting_inds::Vector{SigIndex},
                                regular_up_to::Integer,
                                timer::Timings) where {Char, N}
    
    @info "checking known syzygies, regular up to $(regular_up_to)"
    found_zd = false
    zd_coeffs = Coeff[]
    zd_mons_hsh = MonIdx[]
    zd_ind = zero(SigIndex)

    to_del = Int[]

    @inbounds for (i, (idx, proc_info)) in enumerate(syz_queue)
        syz_mask = basis.syz_masks[idx]
        syz_mon = basis.syz_sigs[idx]
        syz_ind = index(syz_mask)
        tr_ind = tr.syz_ind_to_mat[idx]

        if all(idx -> get(proc_info, idx, false), splitting_inds)
            push!(to_del, i)
        end

        for cofac_ind in reverse(splitting_inds)
            get(proc_info, cofac_ind, false) && continue
            if ind_order.ord[cofac_ind] <= regular_up_to
                proc_info[cofac_ind] = true
                continue
            end
            tim = @elapsed cofac_coeffs, cofac_mons_hsh = construct_module_wrap((syz_ind, syz_mon), basis,
                                                                                basis_ht,
                                                                                tr_ind,
                                                                                tr, char,
                                                                                ind_order, cofac_ind)
            timer.module_time += tim
            proc_info[cofac_ind] = true
            if isempty(cofac_coeffs)
                continue
            end
            if any(gb -> !iszero(my_normal_form(cofac_mons_hsh, cofac_coeffs, basis_ht, gb)), lc_set.gbs)
                found_zd = true
                zd_coeffs, zd_mons_hsh = cofac_coeffs, cofac_mons_hsh
                zd_ind = cofac_ind
                break
            else
                continue
            end
        end
        if found_zd
            break
        end    
    end

    deleteat!(syz_queue, to_del)

    return found_zd, zd_coeffs, zd_mons_hsh, zd_ind
end

# ------------------------ #
# --- helper functions --- #
# ------------------------ #

# write stuff from index i in basis1 to index j in basis2
function overwrite!(basis1::Basis,
                    basis2::Basis,
                    i::Int, j::Int)

    @inbounds begin
        basis2.sigs[j]          = basis1.sigs[i]
        basis2.sigmasks[j]      = basis1.sigmasks[i]
        basis2.sigratios[j]     = basis1.sigratios[i]

        rnodes = copy(basis1.rewrite_nodes[i+1])
        basis2.rewrite_nodes[j+1] = rnodes

        basis2.lm_masks[j]      = basis1.lm_masks[i]

        # TODO: this does not copy, is that safe?
        basis2.monomials[j]     = basis1.monomials[i]
        basis2.coefficients[j]  = basis1.coefficients[i]
        basis2.is_red[j]        = basis1.is_red[i]
    end
end

# make saturators incompatible
function make_sat_incompat!(tags::Tags, ind_order::IndOrder)
    for idx in keys(ind_order.ord)
        gettag(tags, idx) != :sat && continue
        for idx2 in keys(ind_order.ord)
            (gettag(tags, idx2) != :sat || idx2 <= idx) && continue
            ind_order.incompat[(idx, idx2)] = true
        end
    end
end

# homogenize w.r.t. the last variable
function homogenize(F::Vector{P}) where {P <: MPolyRingElem}
    R = parent(first(F))
    S, vars = polynomial_ring(base_ring(R), ["x$i" for i in 1:nvars(R)+1])
    res = typeof(first(F))[]
    for f in F
        ctx = MPolyBuildCtx(S)
        d = total_degree(f)
        for (e, c) in zip(exponent_vectors(f), coefficients(f))
            enew = push!(e, d - sum(e))
            push_term!(ctx, c, e)
        end
        push!(res, finish(ctx))
    end
    return res
end

function is_homog(f)
    d = total_degree(f)
    return all(e -> sum(e) == d, exponent_vectors(f))
end

function new_ind_order(basis::Basis)
    return IndOrder((SigIndex).(collect(1:basis.input_load)),
                    Dict{Tuple{SigIndex, SigIndex}, Bool}(),
                    SigIndex(basis.input_load))
end


#---------------- for debugging --------------------#
function print_sequence(basis::Basis{N},
                        basis_ht::MonomialHashtable,
                        ind_order::IndOrder,
                        lc_set::LocClosedSet,
                        tags::Tags) where N

    inds = sort(collect(1:basis.input_load), by = i -> ind_order.ord[i])
    for i in inds
        mns_hsh = basis.monomials[i]
        mns = [basis_ht.exponents[m] for m in mns_hsh]
        cfs = basis.coefficients[i]
        sig = basis.sigs[i]
        lm = convert_to_pol(ring(lc_set), [first(mns)], [one(Coeff)])
        println("$(i) ------> $(index(basis.sigs[i])), $(lm), $(gettag(tags, index(basis.sigs[i])))")
    end
    for g in lc_set.ineqns
        println(Nemo.leading_monomial(g))
    end
    println("----")
end
        
# test against msolve
function _is_gb(gb::Vector{P}) where {P <: MPolyRingElem}
    gb_pols = gb
    gb_msolve = groebner_basis(Ideal(gb_pols), complete_reduction = true)
    
    lms_gb = (Nemo.leading_monomial).(gb_pols)
    lms_msolve = (Nemo.leading_monomial).(gb_msolve)
    res1 = all(u -> any(v -> divides(u, v)[1], lms_gb), lms_msolve)
    res2 = all(u -> any(v -> divides(u, v)[1], lms_msolve), lms_gb)
    if !(res1 && res2)
        idx = findfirst(u -> all(v -> !divides(u, v)[1], lms_gb), lms_msolve)
        println(lms_gb[idx])
        println(lms_msolve)
        return false
    end
    return true
end

function _is_gb(gb::Vector{Tuple{Tuple{Int, P}, P}}) where {P <: MPolyRingElem}
    gb_pols = [p[2] for p in gb]
    return _is_gb(gb_pols)
end

function _is_gb(basis::Basis{N}, basis_ht::MonomialHashtable,
                tags::Tags, char::Val{Char}) where {N, Char}
    
    R, vrs = polynomial_ring(GF(Int(Char)), ["x$i" for i in 1:N],
                             ordering = :degrevlex)
    G0 = [((Int(index(basis.sigs[i])), convert_to_pol(R, [monomial(basis.sigs[i])], [one(Coeff)])),
           convert_to_pol(R, [basis_ht.exponents[midx]
                             for midx in basis.monomials[i]],
                         basis.coefficients[i]))
          for i in basis.basis_offset:basis.basis_load
          if gettag(tags, index(basis.sigs[i])) != :col]
    seq = [convert_to_pol(R, [basis_ht.exponents[midx] for midx in basis.monomials[i]],
                          basis.coefficients[i])
           for i in 1:basis.input_load
           if gettag(tags, i) != :col]
    if !_is_gb(G0)
        f = open("/tmp/bad.txt", "w+")
        println(f, seq)
        close(f)
        error("no GB")
    end
    return true
end

function new_timer()
    return Timings(0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0)
end

function Base.show(io::IO, timer::Timings)
    @printf io "\n"
    @printf io "symbolic pp:         %.2f\n" timer.sym_pp_time
    @printf io "linear algebra:      %.2f\n" timer.lin_alg_time
    @printf io "select:              %.2f\n" timer.select_time
    @printf io "update:              %.2f\n" timer.update_time
    !iszero(timer.module_time) && @printf io "module construction: %.2f\n" timer.module_time
    !iszero(timer.comp_lc_time) && @printf io "splitting:           %.2f\n" timer.comp_lc_time
end
