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
include("module.jl")
include("normalform.jl")
include("affine_cells.jl")
include("interfaces.jl")


#---------------- user functions --------------------#

@doc Markdown.doc"""
    sig_groebner_basis(sys::Vector{T}; info_level::Int = 0; degbound::Int = 0) where {T <: MPolyRingElem}

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
function sig_groebner_basis(sys::Vector{T}; info_level::Int=0, degbound::Int=0) where {T <: MPolyRingElem}

    # data structure setup/conversion
    sys_mons, sys_coeffs, basis_ht, char, shift = input_setup(sys)

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
        siggb!(basis, pairset, basis_ht, char, shift,
               tags, ind_order, tr, degbound=degbound)
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
        for idx2 in idx+1:length(full_sys)
            ind_order.incompat[(idx, idx2)] = true
        end
    end

    # compute divmasks
    fill_divmask!(basis_ht)
    @inbounds for i in 1:sysl
        basis.lm_masks[i] = basis_ht.hashdata[basis.monomials[i][1]].divmask
    end

    logger = ConsoleLogger(stdout, info_level == 0 ? Warn : Info)
    with_logger(logger) do
        added_unit = siggb!(basis, pairset, basis_ht, char, shift,
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
        timer = Timings(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        lc_sets = sig_decomp!(basis, pairset, basis_ht, char, shift,
                              tags, ind_order, tr, R, timer)
        @info timer
        return lc_sets
    end
end

function sig_kalk_decomp(sys::Vector{T}; info_level::Int=0) where {T <: MPolyRingElem}

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
        timer = Timings(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        lc_sets = sig_kalk_decomp!(basis, pairset, basis_ht, char, shift, tags, R, timer)
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
                tr::Tracer;
                trace::Bool=false,
                degbound::Int=0) where {N, Char, Shift}

    # fake syz queue
    syz_queue = Int[]

    sort_pairset_by_degree!(pairset, 1, pairset.load-1)

    while !iszero(pairset.load)
        if !iszero(degbound) && first(pairset.elems).deg > degbound
            break
        end
	matrix = initialize_matrix(Val(N))
        symbol_ht = initialize_secondary_hash_table(basis_ht)

        _, compat_ind = select_normal!(pairset, basis, matrix,
                                       basis_ht, symbol_ht, ind_order, tags)
        symbolic_pp!(basis, matrix, basis_ht, symbol_ht,
                     ind_order, tags, compat_ind)
        finalize_matrix!(matrix, symbol_ht, ind_order)
        iszero(matrix.nrows) && continue
        tr_mat = echelonize!(matrix, tags, ind_order, char,
                             shift, trace = trace)

        trace && push!(tr.mats, tr_mat)

        added_unit = update_siggb!(basis, matrix, pairset, symbol_ht,
                                   basis_ht, ind_order, tags,
                                   tr, char, syz_queue; trace = trace)
        added_unit && return true
        sort_pairset_by_degree!(pairset, 1, pairset.load-1)
    end

    return false
end

#---------------- functions for splitting --------------------#

function sig_decomp!(basis::Basis{N},
                     pairset::Pairset,
                     basis_ht::MonomialHashtable,
                     char::Val{Char},
                     shift::Val{Shift},
                     tags::Tags,
                     ind_order::IndOrder,
                     tr::Tracer,
                     R::MPolyRing,
                     timer::Timings) where {N, Char, Shift}

    # compute ideal
    eqns = [convert_to_pol(R, [basis_ht.exponents[mdx] for mdx in basis.monomials[i]],
                           basis.coefficients[i])
            for i in 1:basis.input_load]
    X = LocClosedSet{eltype(eqns)}(eqns, [last(gens(R))])
    
    queue = [(basis, pairset, tags, ind_order, X, Int(basis.input_load), Int[], tr)]
    result = LocClosedSet[]

    while !isempty(queue)
        bs, ps, tgs, ind_ord, lc_set, c, syz_queue, tr = popfirst!(queue)
        neqns = num_eqns(lc_set)
        @info "starting component, $(length(queue)) remaining, $(neqns) equations"
        if is_empty_set(lc_set)
            @info "empty component"
            @info "------------------------------------------"
            continue
        elseif codim(lc_set) > c
            @info "superflous component"
            @info "------------------------------------------"
            continue
        elseif codim(lc_set) == neqns
            deleteat!(lc_set.eqns, findall(lc_set.eqns_is_red))
            deleteat!(lc_set.eqns_is_red, findall(lc_set.eqns_is_red)) 
            @info "finished component codim $c"
            push!(result, lc_set)
            @info "------------------------------------------"
            continue
        end
        found_zd, isempt, zd_coeffs,
        zd_mons, zd_ind, _ = siggb_for_split!(bs, ps,
                                              tgs, ind_ord,
                                              basis_ht, tr,
                                              syz_queue,
                                              char, shift, [lc_set],
                                              timer,
                                              maintain_nf=false)
        if found_zd
            @info "splitting component"
            tim = @elapsed bs2, ps2, tgs2,
                           ind_ord2, lc_set2, tr2 = split!(bs, basis_ht,
                                                           zd_mons, zd_coeffs,
                                                           tr, ps,
                                                           zd_ind, tgs,
                                                           ind_ord,
                                                           lc_set)
            timer.comp_lc_time += tim
            pushfirst!(queue, (bs2, ps2, tgs2, ind_ord2, lc_set2, min(c, neqns-1), Int[], tr2))
            pushfirst!(queue, (bs, ps, tgs, ind_ord, lc_set, c, syz_queue, tr))
        else
            deleteat!(lc_set.eqns, findall(lc_set.eqns_is_red))
            deleteat!(lc_set.eqns_is_red, findall(lc_set.eqns_is_red))
            @info "finished component"
            push!(result, lc_set)
        end
        @info "------------------------------------------"
    end
    return result
end

function sig_kalk_decomp!(basis::Basis{N},
                          pairset::Pairset,
                          basis_ht::MonomialHashtable,
                          char::Val{Char},
                          shift::Val{Shift},
                          tags::Tags,
                          R::MPolyRing,
                          timer::Timings) where {N, Char, Shift}

    # index order
    ind_order = new_ind_order(basis)

    # compute ideal
    eqns = [convert_to_pol(R, [basis_ht.exponents[mdx] for mdx in basis.monomials[i]],
                           basis.coefficients[i])
            for i in 1:basis.input_load]
    X = LocClosedSet{eltype(eqns)}(eqns, eltype(eqns)[])
    tracer = new_tracer()
    
    queue = [(basis, pairset, tags, ind_order, [X], Int[], tracer)]
    result = LocClosedSet[]

    while !isempty(queue)
        bs, ps, tgs, ind_ord, lc_sets, syz_queue, tr = popfirst!(queue)
        isempty(lc_sets) && continue
        neqns = num_eqns(first(lc_sets))
        @info "starting set of components, $(length(queue)) remaining, $(neqns) equations, $(length(lc_sets)) components"
        @info "checking for empty components"
        filter!(!is_empty_set, lc_sets)
        if isempty(lc_sets)
            @info "all components empty"
            @info "------------------------------------------"
            continue
        end
        @info "checking for equidimensional components"
        equidim_inds = findall(X -> codim(X) == num_eqns(X), lc_sets)
        @info "$(length(equidim_inds)) components already equidimensional"
        append!(result, lc_sets[equidim_inds])
        # what syntax!
        if !isempty(equidim_inds)
            lc_sets = lc_sets[findall(i -> !(i in equidim_inds), 1:length(lc_sets))]
            if isempty(lc_sets)
                @info "------------------------------------------"
                continue
            end
        end
        
        found_zd, isemptyset, zd_coeffs,
        zd_mons, zd_ind, nz_nf_inds = siggb_for_split!(bs, ps,
                                                       tgs, ind_ord,
                                                       basis_ht, tr, syz_queue,
                                                       char, shift, lc_sets,
                                                       timer, maintain_nf = false)
        if found_zd
            @info "splitting components"
            tim = @elapsed new_lc_sets1, bs2, ps2, tgs2,
                           ind_ord2, new_lc_sets2, tr2 = kalksplit!(bs, basis_ht,
                                                                    zd_mons, zd_coeffs,
                                                                    zd_ind, nz_nf_inds,
                                                                    tags, ind_ord,
                                                                    lc_sets)
            timer.comp_lc_time += tim
            lc_sets1 = vcat(lc_sets[findall(i -> !(i in nz_nf_inds), 1:length(lc_sets))], new_lc_sets1)
            lc_sets2 = new_lc_sets2
            pushfirst!(queue, (bs2, ps2, tgs2, ind_ord2, lc_sets2, Int[], tr2))
            pushfirst!(queue, (bs, ps, tgs, ind_ord, lc_sets1, syz_queue, tr))
        else
            [deleteat!(X.eqns, findall(X.eqns_is_red)) for X in lc_sets]
            @info "finished components $(neqns) eqns"
            append!(result, lc_sets)
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
                          tr::Tracer,
                          syz_queue::Vector{Int},
                          char::Val{Char},
                          shift::Val{Shift},
                          lc_sets::Vector{LocClosedSet{T}},
                          timer::Timings;
                          maintain_nf::Bool=false) where {N, Char, Shift, T <: MPolyRingElem}

    sort_pairset_by_degree!(pairset, 1, pairset.load-1)

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
        tim = @elapsed tr_mat = echelonize!(matrix, tags, ind_order, char, shift)
        timer.lin_alg_time += tim

        push!(tr.mats, tr_mat)

        time = @elapsed update_siggb!(basis, matrix, pairset,
                                      symbol_ht, basis_ht,
                                      ind_order, tags,
                                      tr, char, syz_queue)
        for X in lc_sets
            X.eqns_is_red = basis.is_red[1:basis.input_load]
        end
        timer.update_time += tim

        # check to see if we can split with one of the syzygies
        filter!(idx -> !basis.is_red[index(basis.syz_masks[idx])], syz_queue)
        sort!(syz_queue, by = sz -> (index(basis.syz_masks[sz]), basis.syz_sigs[sz]), lt = lt_squeue)
        if !isempty(syz_queue)
            frst_syz_ind = index(basis.syz_masks[first(syz_queue)])
            poss_syz_new_deg = deg - basis.degs[frst_syz_ind]
            if basis.syz_sigs[first(syz_queue)].deg <= poss_syz_new_deg
                does_split, cofac_coeffs, cofac_mons,
                cofac_ind, nz_nf_inds = process_syz_for_split!(syz_queue, basis_ht,
                                                               basis, tr, ind_order, char, lc_sets,
                                                               tags, timer, 
                                                               maintain_nf = maintain_nf)
            else
                does_split = false
            end
        else
            does_split = false
        end

        if does_split
            return true, false, cofac_coeffs,
                   cofac_mons, cofac_ind, nz_nf_inds
        end

        sort_pairset_by_degree!(pairset, 1, pairset.load-1)
    end
    does_split, cofac_coeffs, cofac_mons,
    cofac_ind, nz_nf_inds = process_syz_for_split!(syz_queue, basis_ht,
                                                   basis, tr, ind_order, char, lc_sets,
                                                   tags, timer, 
                                                   maintain_nf = maintain_nf)

    if does_split
        return true, false, cofac_coeffs, cofac_mons,
               cofac_ind, nz_nf_inds
    end

    return false, false, Coeff[], Monomial{N}[], zero(SigIndex), Int[]
end

function split!(basis::Basis{N},
                basis_ht::MonomialHashtable{N},
                cofac_mons_hsh::Vector{MonIdx},
                cofac_coeffs::Vector{Coeff},
                tr::Tracer,
                pairset::Pairset{N},
                zd_ind::SigIndex,
                tags::Tags,
                ind_order::IndOrder,
                lc_set::LocClosedSet) where N


    @inbounds begin
        # new polynomial
        cofac_mons = [basis_ht.exponents[midx] for midx in cofac_mons_hsh]
        h = convert_to_pol(ring(lc_set), cofac_mons, cofac_coeffs)

        # 2nd component input data
        sorted_inds = collect(1:basis.input_load)
        to_del = findall(idx -> basis.is_red[idx], sorted_inds)
        push!(to_del, zd_ind)
        unique!(sort!(to_del))
        deleteat!(sorted_inds, to_del)
        sort!(sorted_inds, by = ind -> ind_order.ord[ind])

        sys2_mons = copy(basis.monomials[sorted_inds])
        sys2_coeffs = copy(basis.coefficients[sorted_inds])
        lc_set2 = deepcopy(lc_set)
        lc_set2.eqns = lc_set2.eqns[sorted_inds]
        lc_set2.eqns_is_red = lc_set2.eqns_is_red[sorted_inds]
        add_inequation!(lc_set2, h)

        # build basis/pairset for second new system
        basis2, ps2, tags2, ind_ord2, tr2 = fill_structs!(sys2_mons, sys2_coeffs,
                                                          basis_ht, def_tg=:split)

        ind_ord2 = new_ind_order(basis2)

        # 1st component
        zd_deg = basis_ht.exponents[first(cofac_mons_hsh)].deg

        ge_deg_inds = filter(i -> basis.degs[i] > zd_deg, 1:basis.input_load)
        if isempty(ge_deg_inds)
            ord_ind = ind_order.max_ind + one(SigIndex)
        else
            ord_ind, _ = findmin(i -> ind_order.ord[i], ge_deg_inds)
        end

        # insert zd in system
        add_new_sequence_element!(basis, basis_ht, tr,
                                  cofac_mons_hsh, cofac_coeffs,
                                  ind_order, ord_ind, pairset,
                                  tags)

        add_equation!(lc_set, h)
    end

    return basis2, ps2, tags2, ind_ord2, lc_set2, tr2
end    

function kalksplit!(basis::Basis{N},
                    basis_ht::MonomialHashtable{N},
                    cofac_mons_hsh::Vector{MonIdx},
                    cofac_coeffs::Vector{Coeff},
                    zd_ind::SigIndex,
                    nz_nf_inds::Vector{Int},
                    tags::Tags,
                    ind_order::IndOrder,
                    lc_sets::Vector{LocClosedSet{T}}) where {N, T <: MPolyRingElem}


    @inbounds begin
        cofac_mons = [basis_ht.exponents[midx] for midx in cofac_mons_hsh]
        h = convert_to_pol(ring(first(lc_sets)), cofac_mons, cofac_coeffs)

        # compute hull(X, h) for relevant X
        @info "taking hull"
        lc_sets_new1 = LocClosedSet{T}[]
        lc_sets_new2 = LocClosedSet{T}[]
        for X in lc_sets[nz_nf_inds]
            hll = hull(X, h)
            append!(lc_sets_new1, hll)
            nz = add_inequation(X, h)
            push!(lc_sets_new2, nz)
            deleteat!(nz.eqns, zd_ind)
            deleteat!(nz.eqns_is_red, zd_ind)
        end
                
        # 2nd kind of component
        sys2_mons = copy(basis.monomials[1:basis.input_load])
        sys2_coeffs = copy(basis.coefficients[1:basis.input_load])
        deleteat!(sys2_mons, zd_ind)
        deleteat!(sys2_coeffs, zd_ind)

        # build basis/pairset for second new system
        basis2, ps2, tags2, ind_ord2, tr2 = fill_structs!(sys2_mons, sys2_coeffs,
                                                          basis_ht, def_tg=:split)
        ind_ord2 = new_ind_order(basis2)
    end

    return lc_sets_new1, basis2, ps2, tags2, ind_ord2, lc_sets_new2, tr2
end    

function process_syz_for_split!(syz_queue::Vector{Int},
                                basis_ht::MonomialHashtable,
                                basis::Basis{N},
                                tr::Tracer,
                                ind_order::IndOrder,
                                char::Val{Char},
                                lc_sets::Vector{LocClosedSet{T}},
                                tags::Tags,
                                timer::Timings;
                                maintain_nf::Bool=false) where {Char, N,
                                                                T <: MPolyRingElem}
    
    @info "checking known syzygies"
    found_zd = false
    zd_coeffs = Coeff[]
    zd_mons_hsh = MonIdx[]
    zd_ind = zero(SigIndex)
    nz_nf_inds = Int[]
    
    to_del = Int[]

    ind_info = [(index(basis.sigs[i]), basis.is_red[i]) for i in 1:basis.input_load]
    filter!(tpl -> !tpl[2], ind_info)
    sorted_inds = (tpl -> tpl[1]).(ind_info)
    sort!(sorted_inds, by = ind -> ind_order.ord[ind])

    if maintain_nf
        @assert isone(length(lc_sets))
        mns = [[basis_ht.exponents[midx] for midx in basis.monomials[i]]
               for i in basis.basis_offset:basis.basis_load]
        cfs = [basis.coefficients[i] for i in basis.basis_offset:basis.basis_load]
        gb_lens, gb_cfs, gb_exps = _convert_to_msolve(mns, cfs)
    else
        gb_lens, gb_cfs, gb_exps = Int32[], Int32[], Int32[]
    end
    
    @inbounds for (i, idx) in enumerate(syz_queue)
        syz_mask = basis.syz_masks[idx]
        syz_mon = basis.syz_sigs[idx]
        syz_ind = index(syz_mask)
        tr_ind = tr.syz_ind_to_mat[idx]

        for cofac_ind in reverse(sorted_inds)
            tim = @elapsed cofac_coeffs, cofac_mons_hsh = construct_module((syz_ind, syz_mon), basis,
                                                                           basis_ht,
                                                                           tr_ind,
                                                                           tr, char,
                                                                           ind_order, cofac_ind,
                                                                           gb_lens, gb_cfs, gb_exps,
                                                                           maintain_nf=maintain_nf)
            timer.module_time += tim
            if isempty(cofac_coeffs)
                continue
            end
            iszs = [my_iszero_normal_form(cofac_mons_hsh, cofac_coeffs,
                                          basis_ht, lc_set.gb) for lc_set in lc_sets]
            if all(iszs)
                continue
            else
                found_zd = true
                zd_coeffs, zd_mons_hsh = cofac_coeffs, cofac_mons_hsh
                zd_ind = cofac_ind
                nz_nf_inds = findall(b -> !b, iszs)
                break
            end
        end
        if found_zd
            break
        else
            push!(to_del, i)
        end    
    end

    deleteat!(syz_queue, to_del)
            
    if found_zd
        sort_poly!((zd_coeffs, zd_mons_hsh),
                   by = midx -> basis_ht.exponents[midx],
                   lt = lt_drl, rev = true)
        # normalize cofac coefficients
        normalize_cfs!(zd_coeffs, char)
    end

    return found_zd, zd_coeffs, zd_mons_hsh, zd_ind, nz_nf_inds
end

#---------------- helper functions --------------------#

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

# homogenize w.r.t. the last variable
function homogenize(F::Vector{P}) where {P <: MPolyRingElem}
    R = parent(first(F))
    S, vars = polynomial_ring(base_ring(R), ["x$i" for i in 1:nvars(R)+1],
                             ordering = :degrevlex)
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

function Base.show(io::IO, timer::Timings)
    @printf io "\n"
    @printf io "symbolic pp:         %.2f\n" timer.sym_pp_time
    @printf io "linear algebra:      %.2f\n" timer.lin_alg_time
    @printf io "select:              %.2f\n" timer.select_time
    @printf io "update:              %.2f\n" timer.update_time
    @printf io "module construction: %.2f\n" timer.module_time
    @printf io "splitting:           %.2f\n" timer.comp_lc_time
end
