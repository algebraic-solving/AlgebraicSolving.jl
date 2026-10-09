# Tracer Methods

function new_tracer()
    mats = SigTracerMatrix[]
    basis_ind_to_mat = Vector{Int}(undef, init_basis_size)
    syz_ind_to_mat = Int[]
    return SigTracer(mats, basis_ind_to_mat,
                     syz_ind_to_mat, 0,
                     init_basis_size, false, 0)
end

is_complete(tr::SigTracer) = tr.is_complete

function new_tr_mat(nrows::Int,
                    tr::SigTracer)

    # when replaying we hand back the recorded matrix the caller is working on
    # instead of appending to what we are reading
    is_complete(tr) && return tr.mats[tr.curr_mat]

    diag = Vector{Coeff}(undef, nrows)
    mat_data = Vector{Vector{Tuple{Int, Coeff}}}(undef, nrows)

    res = SigTracerMatrix(Vector{Tuple{Sig, Int, Bool}}(undef, nrows),
                          Dict{Sig, Int}(),
                          Dict{Int, Int}(),
                          Dict{Int, Sig}(),
                          diag,
                          mat_data,
                          zero(Exp),
                          Int[],)
    push!(tr.mats, res)
    return res
end

function add_row!(tr_mat::SigTracerMatrix,
                  sig::Sig,
                  row_ind::Int,
                  parent_ind::Int,
                  is_pivot::Bool)

    tr_mat.rows[row_ind] = (sig, parent_ind, is_pivot)
    tr_mat.sig_to_row[sig] = row_ind

    # allocate a row for the tracer matrix
    # at most we subtract (i-1) other rows
    row_ops = Tuple{Int, Coeff}[]
    tr_mat.col_inds_and_coeffs[row_ind] = row_ops
    tr_mat.row_ind_to_sig[row_ind] = sig
    tr_mat.diagonal[row_ind] = one(Coeff)
end

function store_row_op!(tr_mat::SigTracerMatrix,
                       row_ind::Int,
                       pividx::Int,
                       a::Cbuf)
    
    push!(tr_mat.col_inds_and_coeffs[row_ind], (pividx, a))
end

function store_inver!(tr_mat::SigTracerMatrix,
                      row_ind::Int,
                      inver::Coeff)

    tr_mat.diagonal[row_ind] = inver
end

function store_basis_elem!(tr::SigTracer,
                           new_sig::Sig,
                           bas_ind::Int,
                           bas_sz::Int)

    is_complete(tr) && return

    if tr.size != bas_sz
        tr.size = bas_sz
        resize!(tr.basis_ind_to_mat, tr.size)
    end
    tr_mat = last(tr.mats)
    @inbounds row_ind = tr_mat.sig_to_row[new_sig]
    @inbounds tr_mat.is_basis_row[row_ind] = bas_ind
    @inbounds tr.basis_ind_to_mat[bas_ind] = length(tr.mats)
end

function store_syz!(tr::SigTracer)

    is_complete(tr) && return
    push!(tr.syz_ind_to_mat, length(tr.mats)) 
end

# Only the structure is copied. The coefficients stay shared with the parent.
function copy_tracer(tr::SigTracer)
    mats = Vector{SigTracerMatrix}(undef, length(tr.mats))
    @inbounds for (i, m) in enumerate(tr.mats)
        mats[i] = SigTracerMatrix(copy(m.rows),
                                  copy(m.sig_to_row),
                                  copy(m.is_basis_row),
                                  copy(m.row_ind_to_sig),
                                  m.diagonal,
                                  m.col_inds_and_coeffs,
                                  m.deg,
                                  copy(m.toadd))
    end
    return SigTracer(mats, copy(tr.basis_ind_to_mat), copy(tr.syz_ind_to_mat),
                     tr.load, tr.size, tr.is_complete, tr.curr_mat)
end

function shift_tracer!(tr::SigTracer, shift::Int,
                       old_offset::Int,
                       basis::Basis)

    # a recorded tracer already carries the shifted indices
    is_complete(tr) && return

    if tr.size != basis.basis_size
        tr.size = basis.basis_size
        resize!(tr.basis_ind_to_mat, tr.size)
    end
    for i in basis.basis_load:-1:basis.basis_offset
        tr.basis_ind_to_mat[i] = tr.basis_ind_to_mat[i-shift]
    end

    for mat in tr.mats
        for i in eachindex(mat.rows)
            v = mat.rows[i]
            if v[2] >= old_offset
                mat.rows[i] = (v[1], v[2] + shift, v[3])
            end
        end
        for (i, v) in pairs(mat.is_basis_row)
            if v >= old_offset
                mat.is_basis_row[i] = v + shift
            end
        end
    end
end                

function reset_tracers!(ts::TracerStore)
    ts.ind = 1
    ts.syz_ind = 1
end

function record_syz_split!(ts::TracerStore, found::Bool, queue_ind::Int,
                           cofac_ind::SigIndex, to_del::Vector{Int})
    push!(ts.syz_splits, SyzSplit(found, queue_ind, cofac_ind, copy(to_del)))
    ts.syz_ind += 1
    return nothing
end

function next_syz_split!(ts::TracerStore)
    res = ts.syz_splits[ts.syz_ind]
    ts.syz_ind += 1
    return res
end

is_replaying(ts::TracerStore) = ts.recorded

replay_component(ts::TracerStore) = ts.tracers[ts.ind], ts.ranges[ts.ind]

finish_replay!(ts::TracerStore) = ts.ind += 1

function record_component!(ts::TracerStore, tr::SigTracer, start::Int)
    push!(ts.tracers, tr)
    push!(ts.ranges, start:length(tr.mats))
    ts.ind += 1
    return nothing
end

function mark_recorded!(ts::TracerStore)
    ts.recorded && return nothing
    for tr in ts.tracers
        tr.is_complete = true
        for m in tr.mats
            fill!(m.diagonal, one(Coeff))
            for i in eachindex(m.col_inds_and_coeffs)
                isassigned(m.col_inds_and_coeffs, i) && empty!(m.col_inds_and_coeffs[i])
            end
        end
    end
    ts.recorded = true
    return nothing
end

function construct_matrix!(tr::SigTracer,
                           basis::Basis{N},
                           symbol_ht::MonomialHashtable{N},
                           ht::MonomialHashtable{N},
                           ind_order::IndOrder,
                           index::Int) where N

    tr_mat = tr.mats[index]
    nrows = length(tr_mat.rows)
    mat = initialize_matrix(Val(N), nrows)
    @inbounds for i in 1:nrows
        mat.pivots[i] = 0
    end
    # rebuild the rows in their recorded memory order
    @inbounds for row_ind in 1:nrows
        sig, basis_index, is_pivot = tr_mat.rows[row_ind]
        mult = divide(monomial(sig), monomial(basis.sigs[basis_index]))
        lead = write_to_matrix_row!(mat, basis, basis_index, symbol_ht, ht, mult,
                                    sig, row_ind)
        # restore the pivot marks symbolic_pp! would have set
        resize_pivots!(mat, symbol_ht)
        is_pivot && (mat.pivots[lead] = row_ind)
    end
    mat.nrows = nrows
    # symbolic_pp! would normally have grown the pivots along with the columns
    resize_pivots!(mat, symbol_ht)
    finalize_matrix!(mat, symbol_ht, ind_order)
    return mat
end

# dummy methods if we don't want to trace
is_complete(tr::NoTracer) = false

function new_tr_mat(nrows::Int,
                    tr::NoTracer)

    return NoTracerMatrix()
end

function add_row!(tr_mat::NoTracerMatrix,
                  sig::Sig,
                  row_ind::Int,
                  parent_ind::Int,
                  is_pivot::Bool)

    return
end

function store_row_op!(tr_mat::NoTracerMatrix,
                       row_ind::Int,
                       pividx::Int,
                       a::Cbuf)

    return
end

function store_inver!(tr_mat::NoTracerMatrix,
                      row_ind::Int,
                      inver::Coeff)

    return
end

function store_basis_elem!(tr::NoTracer,
                           new_sig::Sig,
                           bas_ind::Int,
                           bas_sz::Int)

    return
end

function store_syz!(tr::NoTracer)

    return
end

copy_tracer(tr::NoTracer) = tr

function shift_tracer!(tr::NoTracer, shift::Int,
                       old_offset::Int,
                       basis::Basis)

    return
end
