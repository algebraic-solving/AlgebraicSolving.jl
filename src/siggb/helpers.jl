# various helper functions

# initialization and memory management

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
    mod_rep_known = Vector{Vector{Bool}}(undef, basis_size)
    mod_reps = Vector{Vector{Polynomial}}(undef, basis_size)
    syz_sigs = Vector{Monomial{N}}(undef, syz_size)
    syz_masks = Vector{MaskSig}(undef, syz_size)
    basis = Basis(sigs, sigmasks, sigratios, rewrite_nodes,
                  lm_masks, monomials, coeffs, is_red,
                  mod_rep_known, mod_reps,
                  syz_sigs, syz_masks, Exp[],
                  input_length,
                  init_basis_size, 0, input_length,
                  input_length + 1, 0,
                  init_syz_size)

    # root node
    basis.rewrite_nodes[1] = [-1, -1]

    return basis
end

function make_room_new_input_el!(basis::Basis,
                                 pairset::Pairset,
                                 tr::Tracer)

    # this whole block just shifts the basis to the right
    # to make room for new input elements
    @inbounds if basis.input_size <= basis.input_load
        basis.input_size *= 2
        
        shift = basis.input_size
        old_offset = basis.basis_offset
        basis.basis_offset += shift
        basis.basis_load += shift
        resize_basis!(basis)

        if tr.load + shift >= tr.size
            tr.size *= 2
            resize!(tr.basis_ind_to_mat, tr.size)
        end

        # adjusts rewrite nodes at the start
        for i in 1:basis.input_load
            if basis.rewrite_nodes[i+1][1] >= 0
                basis.rewrite_nodes[i+1][3] += shift
            end
        end

        for i in basis.basis_load:-1:basis.basis_offset
            basis.sigs[i] = basis.sigs[i-shift]
            basis.sigmasks[i] = basis.sigmasks[i-shift]
            basis.sigratios[i] = basis.sigratios[i-shift]
            basis.lm_masks[i] = basis.lm_masks[i-shift]
            basis.monomials[i] = basis.monomials[i-shift]
            basis.coefficients[i] = basis.coefficients[i-shift]
            basis.is_red[i] = basis.is_red[i-shift]
            basis.mod_rep_known[i] = basis.mod_rep_known[i-shift]
            basis.mod_reps[i] = basis.mod_reps[i-shift]

            # adjust rewrite tree
            rnodes = basis.rewrite_nodes[i-shift+1]
            if rnodes[2] >= old_offset + 1
                rnodes[2] += shift
            end
            for j in 3:3+rnodes[1]
                rnodes[j] += shift
            end
            basis.rewrite_nodes[i+1] = rnodes
        end

        # adjust tracer
        shift_tracer!(tr, shift, old_offset, basis)

        # adjust pairset
        for i in 1:pairset.load
            p = pairset.elems[i]
            if p.top_index >= old_offset
                pairset.elems[i].top_index += shift
            end
            if p.bot_index >= old_offset
                pairset.elems[i].bot_index += shift
            end
        end
    end
end

@inline function compute_shift(i::Int, to_del::Vector{Int})
    idx = findlast(j -> j <= i, to_del)
    if isnothing(idx)
        return 0
    elseif del_indices[idx] == i
        return -1
    else
        return idx
    end
end

# TODO: figure out how to maintain the rewrite tree
function garbage_collect!(basis::Basis{N},
                          pairset::Pairset{N},
                          tr::Tracer,
                          del_indices::Vector{Int}) where N

    isempty(del_indices) && return

    @inbounds for (k, i) in enumerate(del_indices)
        shift = k
        upb = k < length(del_indices) ? del_indices[k+1]-1 : basis.basis_load
        for j in i+1:upb
            basis.sigs[j] = basis.sigs[j+shift]
            basis.sigmasks[j] = basis.sigmasks[j+shift]
            basis.sigratios[j] = basis.sigratios[j+shift]
            basis.lm_masks[j] = basis.lm_masks[j+shift]
            basis.monomials[j] = basis.monomials[j+shift]
            basis.coefficients[j] = basis.coefficients[j+shift]
            basis.is_red[j] = basis.is_red[j+shift]
            basis.mod_rep_known[j] = basis.mod_rep_known[j+shift]
            basis.mod_reps[j] = basis.mod_reps[j+shift]

            # adjust rewrite tree, THIS BREAKS IT
            basis.rewrite_nodes[j+1] = basis.rewrite_nodes[j+shift+1]

            # adjust tracer
            if typeof(Tracer) == SigTracer
                tracer.basis_ind_to_mat[j] = tracer.basis_ind_to_mat[j+shift]
            end
        end
    end

    @inbounds for i in 1:pairset.load
        p = pairset.elems[i]
        shc = compute_shift(p.top_index, del_indices)
        @assert shc != -1
        p.top_index = i - shc
        shc = compute_shift(p.bot_index, del_indices)
        @assert shc != -1
        p.bot_index = i - shc
    end
end

function resize_basis!(basis::Basis)
    if basis.basis_load >= basis.basis_size
        basis.basis_size *= 2
        resize!(basis.sigs, basis.basis_size)
        resize!(basis.sigmasks, basis.basis_size)
        resize!(basis.sigratios, basis.basis_size)
        resize!(basis.rewrite_nodes, basis.basis_size)
        resize!(basis.lm_masks, basis.basis_size)
        resize!(basis.monomials, basis.basis_size)
        resize!(basis.coefficients, basis.basis_size)
        resize!(basis.is_red, basis.basis_size)
        resize!(basis.mod_rep_known, basis.basis_size)
        resize!(basis.mod_reps, basis.basis_size)
    end
end

function resize_syz!(basis::Basis)
    if basis.syz_load >= basis.syz_size
        basis.syz_size *= 2
        resize!(basis.syz_sigs, basis.syz_size)
        resize!(basis.syz_masks, basis.syz_size)
    end
end

function init_pairset(::Val{N}) where N
    ps = Pairset(Vector{SPair{N}}(undef, init_pair_size),
                 0, init_pair_size)
    return ps
end

function resize_pairset!(pairset::Pairset, nnew::Int)
    if pairset.load + nnew >= pairset.size
          resize!(pairset.elems, max(2 * pairset.size,
                                     pairset.load - nnew))
          pairset.size *= 2
    end
end

function initialize_matrix(::Val{N}) where {N}
    rows = Vector{Vector{MonIdx}}(undef, 0)
    pivots = Vector{Int}(undef, 0)
    pivot_size = 0
    sigs = Vector{Sig{N}}(undef, 0)
    parent_inds = Vector{Int}(undef, 0)
    sig_order = Vector{Int}(undef, 0)
    col2hash = Vector{ColIdx}(undef, 0)
    coeffs = Vector{Vector{Coeff}}(undef, 0)
    toadd = Vector{Int}(undef, 0)

    size = 0
    nrows = 0
    ncols = 0
    toadd_length = 0

    return MacaulayMatrix(rows, pivots, pivot_size,
                          sigs, parent_inds, sig_order,
                          col2hash, coeffs, size, nrows,
                          ncols, toadd, toadd_length)
end
    
# Refresh and initialize matrix for `npairs` elements
function reinitialize_matrix!(matrix::MacaulayMatrix, npairs::Int)
    matrix.size = 2 * npairs
    matrix.pivot_size = 2 * npairs
    resize!(matrix.rows, matrix.size)
    resize!(matrix.pivots, matrix.pivot_size)
    for i in 1:matrix.pivot_size
        matrix.pivots[i] = 0
    end
    resize!(matrix.sigs, matrix.size)
    resize!(matrix.parent_inds, matrix.size)
    resize!(matrix.coeffs, matrix.size)
    resize!(matrix.toadd, matrix.size)
    for i in 1:npairs
        matrix.toadd[i] = 0
    end
    return matrix
end

# resize pivots array if needed
@inline function resize_pivots!(matrix::MacaulayMatrix,
                                symbol_ht::MonomialHashtable)
    if matrix.pivot_size < symbol_ht.load 
        pv_size = matrix.pivot_size
        new_pv_size = 2 * (symbol_ht.load)
        resize!(matrix.pivots, new_pv_size)
        @inbounds for j in pv_size+1:new_pv_size 
            matrix.pivots[j] = 0
        end
        matrix.pivot_size = new_pv_size
    end
end

# misc. small functions

function gettag(tags::Tags, i::Integer)
    return get(tags, SigIndex(i), :seq)
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

function new_ind_order(basis::Basis)
    return IndOrder((SigIndex).(collect(1:basis.input_load)),
                    Dict{Tuple{SigIndex, SigIndex}, Bool}(),
                    SigIndex(basis.input_load))
end

# insert new index
function ins_index!(ind_order::IndOrder,
                    new_ord_ind::Integer)

    @inbounds for i in eachindex(ind_order.ord)
        if ind_order.ord[i] >= new_ord_ind
            ind_order.ord[i] += one(SigIndex)
        end
    end
    push!(ind_order.ord, SigIndex(new_ord_ind))
    ind_order.max_ind += 1
end

@inline function comp_sigratio(basis::Basis, ind1::Int, ind2::Int)

    rat1 = basis.sigratios[ind1]
    rat2 = basis.sigratios[ind2]
    if rat1 == rat2
        return true
    end
    return lt_drl(rat1, rat2)
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

# remove pairs that are rewriteable
function remove_red_pairs!(pairset::Pairset)
    iszero(pairset.load) && return
    j = 0 
    @inbounds for i in 1:pairset.load
        iszero(pairset.elems[i].top_index) && continue
        j += 1
        pairset.elems[j] = pairset.elems[i]
    end
    pairset.load = j 
end

function sort_pairset!(pairset::Pairset, from::Int, sz::Int,
                       mord::ModOrd, ind_ord::IndOrder)
    ordr = if mord == :DPOT
        Base.Sort.ord(isless, p -> p.deg, false, Base.Sort.Forward)
    elseif mord == :POT
        Base.Sort.ord(isless, p -> (ind_ord.ord[index(p.top_sig)],
                                    monomial(p.top_sig).deg),
                      false, Base.Sort.Forward)
    end
    sort!(pairset.elems, from, from+sz, def_sort_alg, ordr) 
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
