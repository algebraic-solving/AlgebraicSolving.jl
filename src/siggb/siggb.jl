# sizes for initialization
const init_ht_size = 2^17
const init_basis_size = 10000
const init_syz_size = 1000
const init_pair_size = 10000
const gc_threshold = 300
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
include("helpers.jl")


#---------------- user functions --------------------#

@doc Markdown.doc"""
function sig_groebner_basis(sys::Vector{T}; info_level::Int=0, degbound::Int=0, mod_ord::Symbol=:POT) where {T <: MPolyRingElem}

Compute a Signature Gröbner basis of the sequence `sys` w.r.t. to the
degree reverse lexicographical monomial ordering and the module order
`mod_ord` underlying the computation. The output is a vector of
`Tuple{Tuple{Int64, T}, T}` where the first element indicates the
signature and the second the underlying polynomial.

**Note**: At the moment only ground fields of characteristic `p`, `p` prime, `p < 2^{31}` are supported.
**Note**: If `mod_ord == :DPOT` then the input generators must be homogeneous.
**Note**: The algorithms behaviour may depend heavily on how the elements in `sys` are sorted.

# Arguments
- `sys::Vector{T} where T <: MpolyElem`: input generators.
- `info_level::Int=0`: info level printout: off (`0`, default), computational details (`1`)
- `degbound::Int=0`: Compute a full Gröbner basis if `0` otherwise only up to degree `degbound`.
- `mod_ord::Symbol=:DPOT`: The module monomial order underlying the computation, either `:POT` (position-over-term, default) or `:DPOT` (degree-position-over-term) .

# Example
```jldoctest
julia> using AlgebraicSolving

julia> R, vars = polynomial_ring(GF(17), ["x$i" for i in 1:4])
(Multivariate polynomial ring in 4 variables over GF(17), FqMPolyRingElem[x1, x2, x3, x4])

julia> F = cyclic(R)
FqMPolyRingElem[x1 + x2 + x3 + x4, x1*x2 + x1*x4 + x2*x3 + x3*x4, x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4, x1*x2*x3*x4 + 16]

julia> Fhom = homogenize(F.gens)
4-element Vector{FqMPolyRingElem}:
 x1 + x2 + x3 + x4
 x1*x2 + x2*x3 + x1*x4 + x3*x4
 x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4
 x1*x2*x3*x4 + 16*x5^4

julia> sig_groebner_basis(Fhom, mod_ord = :DPOT)
7-element Vector{Tuple{Tuple{Int64, FqMPolyRingElem}, FqMPolyRingElem}}:
 ((1, 1), x1 + x2 + x3 + x4)
 ((2, 1), x2^2 + 2*x2*x4 + x4^2)
 ((3, 1), x2*x3^2 + x3^2*x4 + 16*x2*x4^2 + 16*x4^3)
 ((4, 1), x2*x3*x4^2 + x3^2*x4^2 + 16*x2*x4^3 + x3*x4^3 + 16*x4^4 + 16*x5^4)
 ((4, x3), x3^3*x4^2 + x3^2*x4^3 + 16*x3*x5^4 + 16*x4*x5^4)
 ((4, x2), x2*x4^4 + x4^5 + 16*x2*x5^4 + 16*x4*x5^4)
 ((4, x2*x3), x3^2*x4^4 + x2*x3*x5^4 + 16*x2*x4*x5^4 + x3*x4*x5^4 + 15*x4^2*x5^4)
```
"""
function sig_groebner_basis(sys::Vector{T}; info_level::Int=0,
                            degbound::Int=0, mod_ord::Symbol=:POT) where {T <: MPolyRingElem}

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
        timer = new_timer()
        _, arit_ops, _ = siggb!(basis, pairset, basis_ht, char, shift,
                                tags, ind_order, tr, timer, degbound,
                                mod_ord)
        @info "$(arit_ops) total submul's"
        @info timer
    end

    # output
    R = parent(first(sys))
    eltp = typeof(first(sys))
    outp = Tuple{Tuple{Int, eltp}, eltp}[]
    @inbounds for i in basis.basis_offset:basis.basis_load
        basis.is_red[i] && continue
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

@doc Markdown.doc"""
function equidimensional_decomposition(I::Ideal{T}, info_level::Int=0) where {T <: MPolyRingElem}

Given a polynomial ideal `I`, return a list of ideals `dec` s.t.
each ideal in `dec` is equidimensional (i.e. has minimal primes
only of one fixed dimension) and s.t. the radical of `I` equals
the intersection of the radicals of the ideals in `dec`.

**Note**: At the moment only ground fields of characteristic `p`, `p` prime, `p < 2^{31}` are supported.

# Arguments
- `I::Ideal{T} where T <: MpolyElem`: input ideal.
- `info_level::Int=0`: info level printout: off (`0`, default), computational details (`1`)

# Example
```jldoctest
julia> using AlgebraicSolving

julia> R, (x, y, z) = polynomial_ring(GF(65521), ["x", "y", "z"])
(Multivariate polynomial ring in 3 variables over GF(65521), Nemo.FqMPolyRingElem[x, y, z])

julia> I = Ideal([x*y - x*z, x*z^2 - x*z, x^2*z - x*z])
Nemo.FqMPolyRingElem[x*y + 65520*x*z, x*z^2 + 65520*x*z, x^2*z + 65520*x*z]

julia> equidimensional_decomposition(I)
3-element Vector{Ideal{Nemo.FqMPolyRingElem}}:
 Nemo.FqMPolyRingElem[x]
 Nemo.FqMPolyRingElem[z, y]
 Nemo.FqMPolyRingElem[z + 65520, y + 65520, x + 65520]
```
"""
function equidimensional_decomposition(I::Ideal{T};
                                       info_level::Int=0) where {T <: MPolyRingElem}

    F = I.gens
    Fhom = homogenize(F)
    sort!(Fhom, by = p -> total_degree(p))
    cells = _sig_decomp(Fhom, info_level = info_level)
    res = typeof(I)[]
    R = parent(I)
    for cell in cells
        for gb in cell.gbs
            push!(res, Ideal(_dehomogenize(gb, R)))
        end
    end
    return res
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
                timer::Timings,
                degbound::Int=0,
                mod_ord::Symbol=:DPOT) where {N, Char, Shift}

    # fake syz queue
    syz_queue = SyzInfo[]
    arit_ops = 0

    nz_conds = Polynomial[]

    sort_pairset!(pairset, 1, pairset.load-1, mod_ord, ind_order)

    while !iszero(pairset.load)
        curr_ind = index(first(pairset.elems).top_sig)
        if !iszero(degbound) && first(pairset.elems).deg > degbound
            break
        end
	matrix = initialize_matrix(Val(N))
        symbol_ht = initialize_secondary_hash_table(basis_ht)

        tim = @elapsed _, compat_ind, sigind = select_normal!(pairset, basis, matrix,
                                                              basis_ht, symbol_ht,
                                                              ind_order, tags,
                                                              mod_ord)
        timer.select_time += tim
        tim = @elapsed symbolic_pp!(basis, matrix, basis_ht, symbol_ht,
                                    ind_order, tags, sigind, compat_ind,
                                    mod_ord)
        timer.sym_pp_time += tim
        finalize_matrix!(matrix, symbol_ht, ind_order)

        if !iszero(matrix.nrows)
            tim = @elapsed arit_ops_new = echelonize!(matrix, tags, ind_order, char,
                                                      shift, tr)
            arit_ops += arit_ops_new
            timer.lin_alg_time += tim

            tim = @elapsed added_unit = update_siggb!(basis, matrix, pairset, symbol_ht,
                                                      basis_ht, ind_order, tags,
                                                      tr, char, syz_queue, mod_ord)
            timer.update_time += tim
            if added_unit
                return true, arit_ops, nz_conds
            end
            sort_pairset!(pairset, 1, pairset.load-1, mod_ord, ind_order)
        end

        p_idx = iszero(pairset.load) ? zero(SigIndex) : index(first(pairset.elems).top_sig)
        if mod_ord == :POT && (iszero(p_idx) || cmp_ind_str(curr_ind, p_idx, ind_order))

            # possible nonzero condition to append in nondeg computation
            if gettag(tags, curr_ind) == :fndegins
                nz_cfs, nz_mons = Coeff[], MonIdx[]
                for (j, syz_msk) in enumerate(basis.syz_masks[1:basis.syz_load])
                    if index(syz_msk) == curr_ind
                        to_add_cfs, to_add_mns = construct_module((curr_ind, basis.syz_sigs[j]),
                                                                  basis,
                                                                  basis_ht,
                                                                  tr.syz_ind_to_mat[j],
                                                                  tr,
                                                                  char,
                                                                  ind_order, curr_ind)
                        mul_by_coeff!(to_add_cfs, rand(one(Coeff):Coeff(Char-1)),
                                      char)
                        nz_cfs, nz_mons = add_pols(nz_cfs, nz_mons,
                                                   to_add_cfs, to_add_mns, char)
                    end
                end
                is_one((nz_cfs, nz_mons), basis_ht) && continue
                sort_poly!((nz_cfs, nz_mons), by = midx -> basis_ht.exponents[midx],
                           lt = lt_drl, rev = true)
                normalize_cfs!(nz_cfs, char)
                push!(nz_conds, (nz_cfs, nz_mons))
            end

            # minimize à la F5c
            min_idx = iszero(p_idx) ? zero(SigIndex) : curr_ind
            minimize!(basis, basis_ht, min_idx, ind_order, tags)
        end
    end
    return false, arit_ops, nz_conds
end

#---------------- functions for splitting --------------------#

function _sig_decomp(sys::Vector{T}; info_level::Int=0) where {T <: MPolyRingElem}

    # data structure setup/conversion
    sys_mons, sys_coeffs, basis_ht, char, shift = input_setup(sys)
    
    # fill basis, pairset, tags
    sysl = length(sys)
    basis, pairset, tags, ind_order, tr = fill_structs!(sys_mons, sys_coeffs,
                                                        basis_ht, def_tg=:split,
                                                        trace=Val(true))

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
    X = loc_closed_set(eqns)
    
    queue = [(basis, pairset, tags, ind_order, X, SyzInfo[], tr)]
    result = LocClosedSet[]

    while !isempty(queue)
        bs, ps, tgs, ind_ord, lc_set, syz_queue, tr = popfirst!(queue)
        neqns = num_eqns(lc_set)
        filter!(gb -> !(one(R) in gb), lc_set.gbs)
        @info "starting component, $(length(queue)) remaining, $(neqns) equations"
        if add_to_output!(result, lc_set)
            @info "done"
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
                          lc_set::LocClosedSet,
                          timer::Timings) where {N, Char, Shift}

    splitting_inds = [index(basis.sigs[i]) for i in 1:basis.input_load]
    filter!(ind -> gettag(tags, ind) == :split, splitting_inds)
    sort!(splitting_inds, by = ind -> ind_order.ord[ind])

    sort_pairset!(pairset, 1, pairset.load-1, :DPOT, ind_order)

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

        tim = @elapsed deg, _, _ = select_normal!(pairset, basis, matrix,
                                                  basis_ht, symbol_ht, ind_order, tags)
        timer.select_time += tim
        tim = @elapsed symbolic_pp!(basis, matrix, basis_ht, symbol_ht,
                                    ind_order, tags)
        timer.sym_pp_time += tim

        finalize_matrix!(matrix, symbol_ht, ind_order)
        iszero(matrix.nrows) && continue
        tim = @elapsed echelonize!(matrix, tags, ind_order, char, shift, tr)
        timer.lin_alg_time += tim

        time = @elapsed update_siggb!(basis, matrix, pairset,
                                      symbol_ht, basis_ht,
                                      ind_order, tags,
                                      tr, char, syz_queue)
        timer.update_time += tim

        # find minimum syzygy index
        sort!(syz_queue, by = sz -> (basis.syz_sigs[sz[1]].deg,
                                     ind_order.ord[index(basis.syz_masks[sz[1]])]))

        # check to see if we can split with one of the syzygies
        if !isempty(syz_queue)
            fs_syz_idx = first(syz_queue)[1]
            min_syz_deg = basis.syz_sigs[fs_syz_idx].deg
            poss_syz_new_deg = deg - basis.degs[basis.input_load]
            if min_syz_deg <= poss_syz_new_deg
                does_split, cofac_coeffs,
                cofac_mons, cofac_ind = process_syz_for_split!(syz_queue, basis_ht,
                                                               basis, tr, ind_order, char, lc_set,
                                                               tags, splitting_inds,
                                                               timer)
                if does_split
                    return true, false, cofac_coeffs,
                    cofac_mons, cofac_ind
                end
            end
        end

        sort_pairset!(pairset, 1, pairset.load-1, :DPOT, ind_order)
    end
    if !isempty(syz_queue)
        sort!(syz_queue, by = sz -> basis.syz_sigs[sz[1]].deg)
        does_split, cofac_coeffs, cofac_mons,
        cofac_ind = process_syz_for_split!(syz_queue, basis_ht,
                                           basis, tr, ind_order, char, lc_set,
                                           tags, splitting_inds,
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
                lc_set::LocClosedSet) where N


    @inbounds begin
        # new polynomial
        cofac_mons = [basis_ht.exponents[midx] for midx in cofac_mons_hsh]
        h = convert_to_pol(ring(lc_set), cofac_mons, cofac_coeffs)
        @info "splitting index $(zd_ind), degree $(total_degree(h))"

        # component with nonzero condition
        sorted_inds = collect(1:basis.input_load)
        deleteat!(sorted_inds, zd_ind)
        sort!(sorted_inds, by = ind -> ind_order.ord[ind])

        sys2_mons = copy(basis.monomials[sorted_inds])
        sys2_coeffs = copy(basis.coefficients[sorted_inds])
        basis2, ps2, tags2, ind_ord2, tr2 = fill_structs!(sys2_mons, sys2_coeffs,
                                                          basis_ht, def_tg=:split,
                                                          trace=Val(true))

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
        lc_set_nz.seq = lc_set.seq[sorted_inds]
        push!(lc_set_hull.seq, h)

        new_codim_ub = min(lc_set.codim_upper_bound, num_eqns(lc_set) - 1)
        lc_set_nz.codim_upper_bound = new_codim_ub
    end

    return lc_set_hull, basis2, ps2, tags2, ind_ord2, lc_set_nz, tr2
end    

function process_syz_for_split!(syz_queue::Vector{SyzInfo},
                                basis_ht::MonomialHashtable,
                                basis::Basis{N},
                                tr::SigTracer,
                                ind_order::IndOrder,
                                char::Val{Char},
                                lc_set::LocClosedSet,
                                tags::Tags,
                                splitting_inds::Vector{SigIndex},
                                timer::Timings) where {Char, N}
    
    @info "checking known syzygies"
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
            if any(gb -> !iszero(normal_form(cofac_mons_hsh, cofac_coeffs, basis_ht, gb)), lc_set.gbs)
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
