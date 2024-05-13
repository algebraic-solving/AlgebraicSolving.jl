# Tracer Methods

function new_tracer()
    mats = SigTracerMatrix[]
    basis_ind_to_mat = Vector{Int}(undef, init_basis_size)
    syz_ind_to_mat = Int[]
    return SigTracer(mats, basis_ind_to_mat,
                     syz_ind_to_mat, 0,
                     init_basis_size)
end

function new_tr_mat(nrows::Int,
                    tr::SigTracer)

    diag = Vector{Coeff}(undef, nrows)
    mat_data = Vector{Vector{Tuple{Int, Coeff}}}(undef, nrows)

    res = SigTracerMatrix(Dict{Sig, Tuple{Int, Int, Bool}}(),
                          Dict{Int, Int}(),
                          Dict{Int, Sig}(),
                          diag,
                          mat_data)
    push!(tr.mats, res)
    return res
end

function add_row!(tr_mat::SigTracerMatrix,
                  sig::Sig,
                  row_ind::Int,
                  parent_ind::Int)

    tr_mat.rows[sig] = (row_ind, parent_ind)

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
                           bas_ind::Int)

    if tr.load >= tr.size
        tr.size *= 2
        resize!(tr.basis_ind_to_mat, tr.size)
    end
    tr_mat = last(tr.mats)
    row_ind, _ = tr_mat.rows[new_sig]
    tr_mat.is_basis_row[row_ind] = bas_ind
    @inbounds tr.basis_ind_to_mat[bas_ind] = length(tr.mats)
    tr.load += 1
end

function store_syz!(tr::SigTracer)

    push!(tr.syz_ind_to_mat, length(tr.mats)) 
end

function shift_tracer!(tr::SigTracer, shift::Int,
                       old_offset::Int,
                       basis::Basis)

    for i in basis.basis_load:-1:basis.basis_offset
        tr.basis_ind_to_mat[i] = tr.basis_ind_to_mat[i-shift]
    end

    for mat in tr.mats
        for i in keys(mat.rows)
            v = mat.rows[i]
            if v[2] >= old_offset
                mat.rows[i] = (v[1], v[2] + shift)
            end
        end
        for (i, v) in pairs(mat.is_basis_row)
            if v >= old_offset
                mat.is_basis_row[i] = v + shift
            end
        end
    end
end                    

# dummy methods if we don't want to trace
function new_tr_mat(nrows::Int,
                    tr::NoTracer)

    return NoTracerMatrix()
end

function add_row!(tr_mat::NoTracerMatrix,
                  sig::Sig,
                  row_ind::Int,
                  parent_ind::Int)

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
                           bas_ind::Int)

    return
end

function store_syz!(tr::NoTracer)

    return
end

function shift_tracer!(tr::NoTracer, shift::Int,
                       old_offset::Int,
                       basis::Basis)

    return
end
