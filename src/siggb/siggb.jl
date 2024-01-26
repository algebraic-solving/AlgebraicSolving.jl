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


#---------------- user functions --------------------#

@doc Markdown.doc"""
    sig_groebner_basis(sys::Vector{T}; info_level::Int = 0, degbound::Int = 0) where {T <: MPolyRingElem}

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
- `degbound::Int=0`: degree bound for Gröbner basis computation, compute a full Gröbner basis if `0` (default) or only up to degree `d`.

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
    basis = fill_basis!(sys_mons, sys_coeffs, basis_ht)
    nv = ngens(parent(first(sys)))
    pairset = init_pairset(Val(nv))
    tags = Tags()
    @inbounds for i in 1:basis.input_load
        add_unit_pair!(basis, pairset, i, basis.degs[i])
        tags[SigIndex(i)] = :seq
    end

    sysl = length(sys)
    # compute divmasks
    fill_divmask!(basis_ht)
    @inbounds for i in 1:sysl
        basis.lm_masks[i] = basis_ht.hashdata[basis.monomials[i][1]].divmask
    end

    logger = ConsoleLogger(stdout, info_level == 0 ? Warn : Info)
    with_logger(logger) do
        siggb!(basis, pairset, basis_ht, char, shift, tags,
               degbound = degbound)
    end

    @assert _is_gb(basis, basis_ht, tags, char)

    # output
    R = parent(first(sys))
    eltp = typeof(first(sys))
    outp = Tuple{Tuple{Int, eltp}, eltp}[]
    @inbounds for i in basis.basis_offset:basis.basis_load
        gettag(tags, index(basis.sigs[i])) == :col && continue
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

function sig_decomp(sys::Vector{T}; info_level::Int=0) where {T <: MPolyRingElem}

    # data structure setup/conversion
    sys_mons, sys_coeffs, basis_ht, char, shift = input_setup(sys)
    
    # fill basis, pairset, tags
    basis = fill_basis!(sys_mons, sys_coeffs, basis_ht)
    nv = ngens(parent(first(sys)))
    pairset = init_pairset(Val(nv))
    tags = Tags()
    @inbounds for i in 1:basis.input_load
        add_unit_pair!(basis, pairset, i, basis.degs[i])
        tags[SigIndex(i)] = :split
    end

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
        lc_sets = sig_decomp!(basis, pairset, basis_ht, char, shift, tags, R, timer)
        @info timer
        return lc_sets
    end
end

function sig_kalk_decomp(sys::Vector{T}; info_level::Int=0) where {T <: MPolyRingElem}

    # data structure setup/conversion
    sys_mons, sys_coeffs, basis_ht, char, shift = input_setup(sys)
    
    # fill basis, pairset, tags
    basis = fill_basis!(sys_mons, sys_coeffs, basis_ht)
    nv = ngens(parent(first(sys)))
    pairset = init_pairset(Val(nv))
    tags = Tags()
    @inbounds for i in 1:basis.input_load
        add_unit_pair!(basis, pairset, i, basis.degs[i])
        tags[SigIndex(i)] = :split
    end

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
                tags::Tags;
                degbound = 0) where {N, Char, Shift}

    # index order
    ind_order = IndOrder((SigIndex).(collect(1:basis.basis_offset-1)),
                         Dict{Tuple{SigIndex, SigIndex}, Bool}(),
                         SigIndex(basis.basis_offset-1))

    # tracer
    tr = new_tracer()

    # fake syz queue
    syz_queue = Int[]

    # sum of degrees of nonzero conditions
    nz_deg = zero(Exp)
    @inbounds for i in 1:basis.input_load
        if gettag(tags, SigIndex(i)) == :col
            nz_deg += basis.degs[i]
        end
    end

    sort_pairset_by_degree!(pairset, 1, pairset.load-1)

    while !iszero(pairset.load)
        if !iszero(degbound) && first(pairset.elems).deg > degbound
            break
        end
	matrix = initialize_matrix(Val(N))
        symbol_ht = initialize_secondary_hash_table(basis_ht)

        deg = select_normal!(pairset, basis, matrix,
                             basis_ht, symbol_ht, ind_order, tags)
        symbolic_pp!(basis, matrix, basis_ht, symbol_ht,
                     ind_order, tags)
        finalize_matrix!(matrix, symbol_ht, ind_order)
        iszero(matrix.nrows) && continue
        tr_mat = echelonize!(matrix, tags, ind_order, char, shift)

        push!(tr.mats, tr_mat)

        update_siggb!(basis, matrix, pairset, symbol_ht,
                      basis_ht, ind_order, tags,
                      tr, char, syz_queue)
        sort_pairset_by_degree!(pairset, 1, pairset.load-1)
    end
end


#---------------- functions for splitting --------------------#

function sig_decomp!(basis::Basis{N},
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
    
    queue = [(basis, pairset, tags, ind_order, X, Int(basis.input_load))]
    result = LocClosedSet[]

    while !isempty(queue)
        bs, ps, tgs, ind_ord, lc_set, c = popfirst!(queue)
        neqns = length(findall(i -> gettag(tgs, i) == :split, 1:bs.input_load))
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
            @info "finished component codim $c"
            push!(result, lc_set)
            @info "------------------------------------------"
            continue
        end
        syz_queue = Int[]
        tr = new_tracer()
        found_zd, isempty, zd_coeffs,
        zd_mons, zd_ind, _, syz_finished = siggb_for_split!(bs, ps,
                                                            tgs, ind_ord,
                                                            basis_ht, tr, syz_queue,
                                                            char, shift, [lc_set],
                                                            timer)
        if found_zd
            @info "splitting component"
            tim = @elapsed bs1, ps1, tgs1, ind_ord1, lc_set1,
                           bs2, ps2, tgs2, ind_ord2, lc_set2 = split!(bs, basis_ht, zd_mons,
                                                                      zd_coeffs, zd_ind, tgs,
                                                                      ind_ord, syz_finished,
                                                                      lc_set)
            timer.comp_lc_time += tim
            pushfirst!(queue, (bs2, ps2, tgs2, ind_ord2, lc_set2, min(c, neqns-1)))
            pushfirst!(queue, (bs1, ps1, tgs1, ind_ord1, lc_set1, c))
        else
            # TODO: this may be dangerous
            to_del_eqns = findall(i -> bs.is_red[i], 1:neqns)
            neqns -= length(to_del_eqns)
            deleteat!(lc_set.eqns, to_del_eqns)
            @info "finished component $(neqns) eqns, codim $(codim(lc_set))"
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
        neqns = length(findall(i -> gettag(tgs, i) == :split, 1:bs.input_load))
        @info "starting set of components, $(length(queue)) remaining, $(neqns) equations, $(length(lc_sets)) components"
        @info "checking for empty components"
        filter!(!is_empty_set, lc_sets)
        if isempty(lc_sets)
            @info "all components empty"
            @info "------------------------------------------"
            continue
        end
        @info "checking for equidimensional components"
        equidim_inds = findall(X -> codim(X) == length(X.eqns), lc_sets)
        @info "$(length(equidim_inds)) components already equidimensional"
        append!(result, lc_sets[equidim_inds])
        # what syntax!
        if !isempty(equidim_inds)
            lc_sets = lc_sets[1:end .∉ equidim_inds]
            if isempty(lc_sets)
                @info "------------------------------------------"
                continue
            end
        end
        
        found_zd, isemptyset, zd_coeffs,
        zd_mons, zd_ind, nz_nf_inds, _ = siggb_for_split!(bs, ps,
                                                          tgs, ind_ord,
                                                          basis_ht, tr, syz_queue,
                                                          char, shift, lc_sets,
                                                          timer)
        if found_zd
            @info "splitting components"
            tim = @elapsed new_lc_sets1, bs2, ps2, tgs2,
                           ind_ord2, new_lc_sets2 = kalksplit!(bs, basis_ht,
                                                               zd_mons, zd_coeffs,
                                                               zd_ind, nz_nf_inds,
                                                               tags, ind_ord,
                                                               lc_sets)
            timer.comp_lc_time += tim
            lc_sets1 = vcat(lc_sets[1:end .∉ nz_nf_inds], new_lc_sets1)
            lc_sets2 = new_lc_sets2
            pushfirst!(queue, (bs2, ps2, tgs2, ind_ord2, lc_sets2, Int[], new_tracer()))
            pushfirst!(queue, (bs, ps, tgs, ind_ord, lc_sets1, syz_queue, tr))
        else
            to_del_eqns = findall(i -> bs.is_red[i], 1:neqns)
            neqns -= length(to_del_eqns)
            [deleteat!(X.eqns, to_del_eqns) for X in lc_sets]
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
                          timer::Timings) where {N, Char, Shift, T <: MPolyRingElem}

    # syz queue
    syz_finished = collect(1:basis.syz_load)

    # module cache
    mod_cache = ModCache{N}()

    sort_pairset_by_degree!(pairset, 1, pairset.load-1)

    while !iszero(pairset.load)
	matrix = initialize_matrix(Val(N))
        symbol_ht = initialize_secondary_hash_table(basis_ht)

        tim = @elapsed deg = select_normal!(pairset, basis, matrix,
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
        timer.update_time += tim

        # check to see if we can split with one of the syzygies
        sort!(syz_queue, by = sz -> monomial(basis.sigs[sz]).deg)
        does_split, cofac_coeffs, cofac_mons,
        cofac_ind, nz_nf_inds = process_syz_for_split!(syz_queue, syz_finished, basis_ht,
                                                       basis, tr, ind_order, char, lc_sets,
                                                       mod_cache, tags, timer)

        if does_split
            return true, false, cofac_coeffs,
                   cofac_mons, cofac_ind, nz_nf_inds,
                   syz_finished
        end

        sort_pairset_by_degree!(pairset, 1, pairset.load-1)
    end
    does_split, cofac_coeffs, cofac_mons,
    cofac_ind, nz_nf_inds = process_syz_for_split!(syz_queue, syz_finished, basis_ht,
                                                   basis, tr, ind_order, char, lc_sets,
                                                   mod_cache, tags, timer)

    if does_split
        return true, false, cofac_coeffs, cofac_mons,
               cofac_ind, nz_nf_inds, syz_finished
    end

    return false, false, Coeff[], Monomial{N}[], zero(SigIndex), Int[], syz_finished
end

function split!(basis::Basis{N},
                basis_ht::MonomialHashtable{N},
                cofac_mons_hsh::Vector{MonIdx},
                cofac_coeffs::Vector{Coeff},
                zd_ind::SigIndex,
                tags::Tags,
                ind_order::IndOrder,
                syz_finished::Vector{Int},
                lc_set::LocClosedSet) where N


    @inbounds begin
        # 1st component
        sys1_mons = copy(basis.monomials[1:basis.input_load])
        sys1_coeffs = copy(basis.coefficients[1:basis.input_load])
        zd_deg = basis_ht.exponents[first(cofac_mons_hsh)].deg

        nz_from = findfirst(i -> gettag(tags, i) == :col, 1:basis.input_load)
        if isnothing(nz_from)
            nz_from = basis.input_load+1
        end
        ins_ind = findfirst(i -> basis.degs[i] > zd_deg, 1:nz_from-1)
        if isnothing(ins_ind)
            ins_ind = nz_from
        end

        # insert zd in system
        insert!(sys1_mons, ins_ind, cofac_mons_hsh)
        insert!(sys1_coeffs, ins_ind, cofac_coeffs)

        # tags for first system
        tags1 = Tags()
        tags1[ins_ind] = :split
        for k in keys(tags)
            if k >= ins_ind
                tags1[k+1] = tags[k]
            elseif k < ins_ind
                tags1[k] = tags[k]
            end
        end

        basis1 = fill_basis!(sys1_mons, sys1_coeffs,
                             basis_ht)
        ps1 = init_pairset(Val(N))
        @inbounds for i in 1:basis1.input_load
            if gettag(tags1, i) == :col
                basis1.is_red[i] = true
            end
            add_unit_pair!(basis1, ps1, i, basis1.degs[i])
        end

        for idx in syz_finished
            syz_sig_mon = basis.syz_sigs[idx]
            syz_msk = basis.syz_masks[idx]
            syz_sig_idx = index(syz_msk)
            new_syz_idx = syz_sig_idx >= ins_ind ? syz_sig_idx+1 : syz_sig_idx
            resize_syz!(basis1)
            l = basis1.syz_load
            basis1.syz_sigs[l+1] = syz_sig_mon
            basis1.syz_masks[l+1] = (new_syz_idx, mask(syz_msk))
            basis1.syz_load += 1
        end

        cofac_mons = [basis_ht.exponents[midx] for midx in cofac_mons_hsh]
        lc_set1 = deepcopy(lc_set)
        h = convert_to_pol(ring(lc_set1), cofac_mons, cofac_coeffs)
        add_equation!(lc_set1, h)

        ind_ord1 = new_ind_order(basis1)
                
        # 2nd component
        sys2_mons = copy(basis.monomials[1:basis.input_load])
        sys2_coeffs = copy(basis.coefficients[1:basis.input_load])
        push!(sys2_mons, cofac_mons_hsh)
        push!(sys2_coeffs, cofac_coeffs)
        deleteat!(sys2_mons, zd_ind)
        deleteat!(sys2_coeffs, zd_ind)

        # new tags
        tags2 = Tags()
        for k in keys(tags)
            if k < nz_from-1
                tags2[k] = :split
            else
                tags2[k] = :col
            end
        end

        # build basis/pairset for second new system
        basis2 = fill_basis!(sys2_mons, sys2_coeffs,
                             basis_ht)

        ps2 = init_pairset(Val(N))
        @inbounds for i in 1:basis2.input_load
            if gettag(tags2, i) == :col
                basis2.is_red[i] = true
            end
            add_unit_pair!(basis2, ps2, i, basis2.degs[i])
        end

        lc_set2 = deepcopy(lc_set)
        add_inequation!(lc_set2, h)
        deleteat!(lc_set2.eqns, zd_ind)

        ind_ord2 = new_ind_order(basis2)
    end

    return basis1, ps1, tags1, ind_ord1, lc_set1,
           basis2, ps2, tags2, ind_ord2, lc_set2
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
        lc_sets_new1 = LocClosedSet{T}[]
        for X in lc_sets[nz_nf_inds]
            append!(lc_sets_new1, hull(X, h))
        end
                
        # 2nd kind of component
        sys2_mons = copy(basis.monomials[1:basis.input_load])
        sys2_coeffs = copy(basis.coefficients[1:basis.input_load])
        deleteat!(sys2_mons, zd_ind)
        deleteat!(sys2_coeffs, zd_ind)

        nz_from = findfirst(i -> gettag(tags, i) == :col, 1:basis.input_load)
        if isnothing(nz_from)
            nz_from = basis.input_load+1
        end
        # new tags
        tags2 = Tags()
        for k in keys(tags)
            if k < nz_from-1
                tags2[k] = :split
            else
                tags2[k] = :col
            end
        end

        # build basis/pairset for second new system
        basis2 = fill_basis!(sys2_mons, sys2_coeffs,
                             basis_ht)

        ps2 = init_pairset(Val(N))
        @inbounds for i in 1:basis2.input_load
            if gettag(tags2, i) == :col
                basis2.is_red[i] = true
            end
            add_unit_pair!(basis2, ps2, i, basis2.degs[i])
        end

        lc_sets_new2 = [add_inequation(X, h) for X in lc_sets[nz_nf_inds]]
        for X in lc_sets_new2
            deleteat!(X.eqns, zd_ind)
        end

        ind_ord2 = new_ind_order(basis2)
    end

    return lc_sets_new1, basis2, ps2, tags2, ind_ord2, lc_sets_new2
end    

function process_syz_for_split!(syz_queue::Vector{Int},
                                syz_finished::Vector{Int},
                                basis_ht::MonomialHashtable,
                                basis::Basis{N},
                                tr::Tracer,
                                ind_order::IndOrder,
                                char::Val{Char},
                                lc_sets::Vector{LocClosedSet{T}},
                                mod_cache::ModCache{N},
                                tags::Tags,
                                timer::Timings) where {Char, N, T <: MPolyRingElem}
    
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
    
    @inbounds for (i, idx) in enumerate(syz_queue)
        syz_mask = basis.syz_masks[idx]
        syz_mon = basis.syz_sigs[idx]
        syz_ind = index(syz_mask)
        tr_ind = tr.syz_ind_to_mat[idx]

        push!(to_del, i)
        for cofac_ind in reverse(sorted_inds)
            tim = @elapsed cofac_coeffs, cofac_mons_hsh = construct_module((syz_ind, syz_mon), basis,
                                                                           basis_ht,
                                                                           tr_ind,
                                                                           tr, char,
                                                                           ind_order, cofac_ind,
                                                                           mod_cache)
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
            push!(syz_finished, idx)
        end
        found_zd && break

    end

    deleteat!(syz_queue, to_del)
            
    if found_zd
        s = sortperm(zd_mons_hsh, by = midx -> basis_ht.exponents[midx],
                     lt = lt_drl, rev = true)
        zd_mons_hsh = zd_mons_hsh[s]
        zd_coeffs = zd_coeffs[s]
        # normalize cofac coefficients
        inver = inv(first(zd_coeffs), char)
        @inbounds for i in eachindex(zd_coeffs)
            if isone(i)
                zd_coeffs[i] = one(Coeff)
                continue
            end
            zd_coeffs[i] = mul(inver, zd_coeffs[i], char)
        end
    end

    return found_zd, zd_coeffs, zd_mons_hsh, zd_ind, nz_nf_inds
end

#---------------- functions for setting up data structures --------------------#

function input_setup(sys::Vector{<:MPolyRingElem})

    if isempty(sys)
        error("Input system is empty.")
    end
    
    R = first(sys).parent
    Rchar = characteristic(R)

    # check if input is ok
    if Rchar > 2^31 || iszero(Rchar)
        error("At the moment we only support finite fields up to prime characteristic < 2^31.")
    end
    sysl = length(sys)
    degs = Vector{Exp}(undef, sysl)
    @inbounds for (i, f) in enumerate(sys)
        deg = total_degree(f)
        if deg > typemax(Exp)
            error("input degrees too large.")
        end
        degs[i] = Exp(deg)
        for m in exponent_vectors(f)
            if sum(m) != deg
                error("input system must be homogeneous.")
            end
        end
    end

    # constants for fast arithmetic
    char = Val(Coeff(Rchar.d))
    shift = Val(maxshift(char))

    # convert to and initialize our data structures
    nv = nvars(R)
    basis_ht = initialize_basis_hash_table(Val(nv))

    sys_mons = Vector{Vector{MonIdx}}(undef, sysl)
    sys_coeffs = Vector{Vector{Coeff}}(undef, sysl)

    # store initial pols in basis and pairset
    @inbounds for i in 1:sysl
        f = sys[i]
        lf = length(f)

        # gather up monomials and coeffs
        exps = collect(exponent_vectors(f))
        cfs = collect(coefficients(f))
        mons = Vector{MonIdx}(undef, lf)
        coeffs = Vector{Coeff}(undef, lf)
        inver = one(Coeff)
        @inbounds for j in 1:lf
            m = monomial(SVector{nv}((Exp).(exps[j])))
            eidx = insert_in_hash_table!(basis_ht, m)
            if isone(j)
                inver = inv(Coeff(lift(ZZ, cfs[1]).d), char)
            end
            cf = isone(j) ? one(Coeff) : mul(inver, Coeff(lift(ZZ, cfs[j]).d), char)
            mons[j] = eidx
            coeffs[j] = cf
        end
        s = sortperm(mons, by = eidx -> basis_ht.exponents[eidx],
                     lt = lt_drl, rev = true)
        mons = mons[s]
        coeffs = coeffs[s]
        sys_mons[i] = copy(mons)
        sys_coeffs[i] = copy(coeffs)
    end

    return sys_mons, sys_coeffs, basis_ht, char, shift
end

function fill_basis!(sys_mons::Vector{Vector{MonIdx}},
                     sys_coeffs::Vector{Vector{Coeff}},
                     basis_ht::MonomialHashtable{N}) where N

    # initialize basis
    sysl = length(sys_mons)
    basis = new_basis(init_basis_size, init_syz_size, sysl, Val(N))

    @inbounds for i in 1:sysl
        s_ind = SigIndex(i)
        mons = sys_mons[i]
        coeffs = sys_coeffs[i]
        lm = basis_ht.exponents[first(mons)]
        lm_mask = divmask(lm, basis_ht.divmap, basis_ht.ndivbits)

        add_input_element!(basis, s_ind,
                           sys_mons[i], sys_coeffs[i],
                           lm_mask, lm)
    end

    return basis
end

function convert_to_pol(R::MPolyRing,
                        exps::Vector{<:Monomial},
                        coeffs::Vector{Coeff})

    ctx = MPolyBuildCtx(R)
    for (e, c) in zip(exps, coeffs)
        push_term!(ctx, base_ring(R)(c), Vector{Int}(e.exps))
    end
    return finish(ctx)
end

function new_basis(basis_size, syz_size,
                   input_length, ::Val{N}) where N

    sigs = Vector{Sig{N}}(undef, basis_size)
    sigmasks = Vector{MaskSig}(undef, basis_size)
    sigratios = Vector{Monomial{N}}(undef, basis_size)
    rewrite_nodes = Vector{Vector{Int}}(undef, basis_size+1)
    lm_masks = Vector{DivMask}(undef, basis_size)
    monomials = Vector{Vector{MonIdx}}(undef, basis_size)
    coeffs = Vector{Vector{Coeff}}(undef, basis_size)
    is_red = Vector{Bool}(undef, basis_size)
    syz_sigs = Vector{Monomial{N}}(undef, syz_size)
    syz_masks = Vector{MaskSig}(undef, syz_size)
    basis = Basis(sigs, sigmasks, sigratios, rewrite_nodes,
                  lm_masks, monomials, coeffs, is_red,
                  syz_sigs, syz_masks, Exp[],
                  input_length,
                  init_basis_size, 0, input_length,
                  input_length + 1, 0,
                  init_syz_size)

    # root node
    basis.rewrite_nodes[1] = [-1, -1]

    return basis
end

function init_pairset(::Val{N}) where N
    ps = Pairset(Vector{SPair{N}}(undef, init_pair_size),
                 0, init_pair_size)
    return ps
end

function add_input_element!(basis::Basis{N},
                            ind::SigIndex,
                            mons::Vector{MonIdx},
                            coeffs::Vector{Coeff},
                            lm_divm::DivMask,
                            lm::Monomial) where N

    @inbounds begin
        one_mon = one_monomial(Monomial{N})

        # signature
        sig = (ind, one_mon)

        l = basis.input_load + 1

        # store stuff in basis
        basis.sigs[l] = sig
        basis.sigmasks[l] = (ind, zero(DivMask))
        basis.sigratios[l] = lm
        basis.rewrite_nodes[l+1] = [-1, 1]
        basis.monomials[l] = mons
        basis.coefficients[l] = coeffs
        basis.is_red[l] = false
        push!(basis.degs, lm.deg)
        basis.lm_masks[l] = lm_divm
        basis.input_load += 1

        # add child to rewrite root
        push!(basis.rewrite_nodes[1], l+1)
        basis.rewrite_nodes[1][1] += 1
    end
end

function add_unit_pair!(basis::Basis{N},
                        pairset::Pairset{N},
                        ind::Integer,
                        deg::Exp) where N

    @inbounds basis.is_red[ind] && return
    one_mon = one_monomial(Monomial{N})
    zero_sig = (zero(SigIndex), one_mon)
    @inbounds sig = basis.sigs[ind]
    @inbounds msk = basis.sigmasks[ind]
    pairset.elems[pairset.load+1] = SPair{N}(sig, zero_sig, mask(msk),
                                             zero(DivMask), Int(ind),
                                             0, deg)
    pairset.load += 1
    return
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

function new_ind_order(basis::Basis)
    return IndOrder((SigIndex).(collect(1:basis.input_load)),
                    Dict{Tuple{SigIndex, SigIndex}, Bool}(),
                    SigIndex(basis.input_load))
end


#---------------- for debugging --------------------#
function print_sequence(basis::Basis{N},
                        basis_ht::MonomialHashtable,
                        ind_order::IndOrder,
                        tags::Tags) where N

    inds = sort(collect(1:basis.input_load), by = i -> ind_order.ord[i])
    for i in inds
        mns_hsh = basis.monomials[i]
        mns = [basis_ht.exponents[m] for m in mns_hsh]
        cfs = basis.coefficients[i]
        sig = basis.sigs[i]
        println("$(first(mns)) ------> $(index(basis.sigs[i])), $(gettag(tags, index(basis.sigs[i]))), $(basis.degs[i])")
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
    return res1 && res2
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
