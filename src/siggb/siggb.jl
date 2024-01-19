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
    tags[basis.input_load] = :col

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
        push_term!(ctx, 1, Vector{Int}(monomial(s).exps))
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
        lc_sets = sig_decomp!(basis, pairset, basis_ht, char, shift, tags, R)
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
                     R::MPolyRing) where {N, Char, Shift}

    # index order
    ind_order = IndOrder((SigIndex).(collect(1:basis.basis_offset-1)),
                         Dict{Tuple{SigIndex, SigIndex}, Bool}(),
                         SigIndex(basis.basis_offset-1))

    # compute ideal
    eqns = [convert_to_pol(R, [basis_ht.exponents[mdx] for mdx in basis.monomials[i]],
                           basis.coefficients[i])
            for i in 1:basis.input_load]
    X = LocClosedSet{eltype(eqns)}(eqns, eltype(eqns)[])
    
    queue = [(basis, pairset, tags, ind_order, X)]
    result = LocClosedSet[]

    while !isempty(queue)
        bs, ps, tgs, ind_ord, lc_set = popfirst!(queue)
        @info "starting component"
        if is_empty_set(lc_set)
            @info "empty/superflous component"
            @info "------------------------------------------"
            continue
        end
        found_zd, isempty, zd_coeffs,
        zd_mons, zd_ind, syz_finished = siggb_for_split!(bs, ps,
                                                         tgs, ind_ord,
                                                         basis_ht, char,
                                                         shift, lc_set)
        if found_zd
            @info "splitting component"
            bs1, ps1, tgs1, ind_ord1, lc_set1,
            bs2, ps2, tgs2, ind_ord2, lc_set2 = split!(bs, basis_ht, zd_mons,
                                                       zd_coeffs, zd_ind, tags,
                                                       ind_ord, lc_set)
            pushfirst!(queue, (bs2, ps2, tgs2, ind_ord2, lc_set2))
            pushfirst!(queue, (bs1, ps1, tgs1, ind_ord1, lc_set1))
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
                          char::Val{Char},
                          shift::Val{Shift},
                          lc_set::LocClosedSet) where {N, Char, Shift}

    # tracer
    tr = new_tracer()

    # syz queue
    syz_queue = Int[]
    syz_finished = collect(1:basis.syz_load)

    sort_pairset_by_degree!(pairset, 1, pairset.load-1)

    while !iszero(pairset.load)
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

        update_siggb!(basis, matrix, pairset,
                      symbol_ht, basis_ht,
                      ind_order, tags,
                      tr, char, syz_queue)

        # check to see if we can split with one of the syzygies
        sort!(syz_queue, by = sz -> monomial(basis.sigs[sz]).deg)
        does_split, cofac_coeffs, cofac_mons,
        cofac_ind = process_syz_for_split!(syz_queue, syz_finished, basis_ht,
                                           basis, tr, ind_order, char, lc_set)

        if does_split
            return true, false, cofac_coeffs,
                   cofac_mons, cofac_ind, syz_finished
        end

        sort_pairset_by_degree!(pairset, 1, pairset.load-1)
    end
    does_split, cofac_coeffs, cofac_mons,
    cofac_ind = process_syz_for_split!(syz_queue, syz_finished, basis_ht,
                                       basis, tr, ind_order, char, lc_set)

    if does_split
        return true, false, cofac_coeffs, cofac_mons,
               cofac_ind, syz_finished
    end

    return false, false, Coeff[], Monomial{N}[], zero(SigIndex), syz_finished
end

function split!(basis::Basis{N},
                basis_ht::MonomialHashtable,
                cofac_mons::Vector{Monomial{N}},
                cofac_coeffs::Vector{Coeff},
                zd_ind::SigIndex,
                tags::Tags,
                ind_order::IndOrder,
                lc_set::LocClosedSet) where N


    @inbounds begin
        # 1st component
        sys1_mons = copy(basis.monomials[1:basis.input_load])
        sys1_coeffs = copy(basis.coefficients[1:basis.input_load])
        
        # find out where to insert zero divisor
        zd_deg = first(cofac_mons).deg
        split_inds = filter(ind -> gettag(tags, SigIndex(ind)) == :split,
                            (index).(basis.sigs[1:basis.input_load]))
        split_inds_deg_geq = filter(ind -> basis.degs[ind] > zd_deg, split_inds)
        if isempty(split_inds_deg_geq)
            ins_ind, _ = findmax(ind -> ind_order.ord[ind], split_inds)
            ins_ind += one(SigIndex)
        else
            ins_ind, _ = findmin(ind -> ind_order.ord[ind], split_inds_deg_geq)
        end
        ind_ord1 = deepcopy(ind_order)
        ins_index!(ind_ord1, ins_ind)

        # insert zd in system
        cofac_mons_hashed = [insert_in_hash_table!(basis_ht, m) for m in cofac_mons]
        push!(sys1_mons, cofac_mons_hashed)
        push!(sys1_coeffs, cofac_coeffs)

        # build basis/pairset/tags for first new system
        basis1 = fill_basis!(sys1_mons, sys1_coeffs,
                             basis_ht)
        tags1 = copy(tags)
        tags1[basis1.input_load] = :split
        ps1 = init_pairset(Val(N))
        @inbounds for i in 1:basis1.input_load
            if i <= basis.input_load && basis.is_red[i]
                basis1.is_red[i] = true
            end
            add_unit_pair!(basis1, ps1, i, basis1.degs[i])
        end

        lc_set1 = deepcopy(lc_set)
        h = convert_to_pol(ring(lc_set1), cofac_mons, cofac_coeffs)
        add_equation!(lc_set1, h)

        # 2nd component
        sys2_mons = copy(basis.monomials[1:basis.input_load])
        sys2_coeffs = copy(basis.coefficients[1:basis.input_load])
        ind_ord2 = deepcopy(ind_order)

        # new tags
        tags2 = copy(tags)

        # build basis/pairset for second new system
        basis2 = fill_basis!(sys2_mons, sys2_coeffs,
                             basis_ht)
        basis2.is_red[zd_ind] = true

        ps2 = init_pairset(Val(N))
        @inbounds for i in 1:basis2.input_load
            if i <= basis.input_load && basis.is_red[i]
                basis2.is_red[i] = true
            end
            add_unit_pair!(basis2, ps2, i, basis2.degs[i])
        end

        lc_set2 = deepcopy(lc_set)
        add_ineqation!(lc_set2, h)
    end

    # print_sequence(basis, basis_ht, ind_order, tags, ring(lc_set))
    # print_sequence(basis1, basis_ht, ind_ord1, tags1, ring(lc_set1))
    # print_sequence(basis2, basis_ht, ind_ord2, tags2, ring(lc_set2))

    return basis1, ps1, tags1, ind_ord1, lc_set1,
           basis2, ps2, tags2, ind_ord2, lc_set2
end    

function process_syz_for_split!(syz_queue::Vector{Int},
                                syz_finished::Vector{Int},
                                basis_ht::MonomialHashtable,
                                basis::Basis{N},
                                tr::Tracer,
                                ind_order::IndOrder,
                                char::Val{Char},
                                lc_set::LocClosedSet) where {Shift, Char, N}
    
    @info "checking known syzygies"
    found_zd = false
    zd_coeffs = Coeff[]
    zd_mons = Monomial{N}[]
    zd_ind = zero(SigIndex)
    
    to_del = Int[]

    sorted_inds = (tpl -> tpl[2]).(sort([(basis.degs[i], index(basis.sigs[i])) for i in 1:basis.input_load],
                                        by = tpl -> tpl[1], rev = true))
    
    @inbounds for (i, idx) in enumerate(syz_queue)
        syz_mask = basis.syz_masks[idx]
        syz_mon = basis.syz_sigs[idx]
        syz_ind = index(syz_mask)
        tr_ind = tr.syz_ind_to_mat[idx]

        for cofac_ind in sorted_inds
            mod_rep = construct_module((syz_ind, syz_mon), basis, tr_ind,
                                       tr, char,
                                       ind_order.max_ind,
                                       ind_order, cofac_ind)
            cofac_mons, cofac_coeffs = mod_rep[cofac_ind][2],
            mod_rep[cofac_ind][1]
            if isempty(cofac_coeffs)
                continue
            end
            isz = my_iszero_normal_form(cofac_mons, cofac_coeffs, lc_set.gb)
            if isz
                continue
            else
                found_zd = true
                zd_coeffs, zd_mons = cofac_coeffs, cofac_mons
                zd_ind = cofac_ind
                break
            end
        end

        push!(syz_finished, idx)
    end
            
    if found_zd
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

    return found_zd, zd_coeffs, zd_mons, zd_ind
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
                inver = inv(Coeff(cfs[1].data), char)
            end
            cf = isone(j) ? one(Coeff) : mul(inver, Coeff(cfs[j].data), char)
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
        push_term!(ctx, c, Vector{Int}(e.exps))
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

#---------------- for locally closet sets --------------------#

# for displaying locally closed sets
function Base.show(io::IO, lc::LocClosedSet)
    string_rep = "V("
    for (i, f) in enumerate(lc.eqns)
        if i != length(lc.eqns)
            string_rep *= "$f, "
        else
            string_rep *= "$(f)) \\ "
        end
    end
    string_rep *= "V("
    for (i, f) in enumerate(lc.ineqns)
        if i != length(lc.ineqns)
            string_rep *= "($f)*"
        else
            string_rep *= "($f)"
        end
    end
    string_rep *= ")"
    print(io, string_rep)
end

function ring(X::LocClosedSet)
    return parent(first(X.eqns))
end

function is_empty_set(X::LocClosedSet)
    R = ring(X)
    return one(R) in X.gb
end

function add_equation!(X::LocClosedSet, f::MPolyRingElem)
    @info "adding equation"
    push!(X.eqns, f)
    X.gb = saturate(push!(X.gb, f), X.ineqns)
end

function add_ineqation!(X::LocClosedSet, h::MPolyRingElem)
    @info "adding inequation"
    push!(X.ineqns, h)
    X.gb = saturate(X.gb, h)
end

function saturate(F::Vector{P}, nzs::Vector{P}) where {P <: MPolyRingElem}
    res = F
    if isempty(nzs)
        res = groebner_basis(Ideal(res), complete_reduction = true)
    else
        for h in nzs
            res = saturate(res, h)
        end
    end
    return res
end

function saturate(F::Vector{P}, nz::P) where {P <: MPolyRingElem}
    R = parent(first(F))
    S, vars = polynomial_ring(base_ring(R), pushfirst!(["x$i" for i in 1:nvars(R)], "t"),
                             ordering = :degrevlex)
    Fconv = typeof(first(F))[]
    for f in F
        ctx = MPolyBuildCtx(S)
        for (e, c) in zip(exponent_vectors(f), coefficients(f))
            enew = pushfirst!(e, 0)
            push_term!(ctx, c, e)
        end
        push!(Fconv, finish(ctx))
    end

    ctx = MPolyBuildCtx(S)
    for (e, c) in zip(exponent_vectors(nz), coefficients(nz))
        enew = pushfirst!(e, 0)
        push_term!(ctx, c, e)
    end
    push!(Fconv, first(vars)*finish(ctx) - 1)

    gb = groebner_basis(Ideal(Fconv), complete_reduction = true, eliminate = 1)

    # convert back to original ring
    res = Vector{P}(undef, length(gb))

    for (i, p) in enumerate(gb)
        ctx = MPolyBuildCtx(R)
        for (e, c) in zip(exponent_vectors(p), coefficients(p))
            push_term!(ctx, c, e[2:end])
        end
        res[i] = finish(ctx)
    end
    return res
end
    
#---------------- for debugging --------------------#
function print_sequence(basis::Basis{N},
                        basis_ht::MonomialHashtable,
                        ind_order::IndOrder,
                        tags::Tags,
                        R::MPolyRing) where N

    inds = sort(collect(1:basis.input_load), by = i -> ind_order.ord[i])
    seq = MPolyRingElem[]
    for i in inds
        basis.is_red[i] && continue
        mns_hsh = basis.monomials[i]
        mns = [basis_ht.exponents[m] for m in mns_hsh]
        cfs = basis.coefficients[i]
        sig = basis.sigs[i]
        p = convert_to_pol(R, mns, cfs)
        push!(seq, p)
        to_pr = length(collect(exponent_vectors(p))) < 5 ? p : Nemo.leading_monomial(p)
        println("$(to_pr) ------> $(index(basis.sigs[i])), $(gettag(tags, index(basis.sigs[i])))")
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
    return _is_gb(G0)
end
