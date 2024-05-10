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
                          tr_diagonal,
                          tr_mat_data)
    push!(tr.mats, res)
    return res
end

function add_row!(tr_mat::SigTracerMatrix,
                  sig::Sig,
                  row_ind::Int,
                  parent_ind::Int,
                  nred::Int)

    tr_mat.rows[row_sig] = (row_ind, parent_ind)

    # allocate a row for the tracer matrix
    # at most we subtract (i-1) other rows
    row_ops = Vector{Tuple{Int, Coeff}}(undef, nred - 1)
    tr_mat.col_inds_and_coeffs[row_ind] = row_ops
    tr_mat.row_ind_to_sig[row_ind] = row_sig
    tr_mat.diagonal[row_ind] = one(Coeff)
end

function resize_tracer_row_ops!(tr_mat::SigTracerMatrix,
                                row_ind::Int,
                                sz::Int)

    resize!(tr_mat.col_inds_and_coeffs[row_ind], sz)
end

function store_row_op!(tr_mat::SigTracerMatrix,
                       row_ind::Int,
                       pividx::Int,
                       a::Coeff)
    
    tr_mat.col_inds_and_coeffs[row_ind][n_row_subs] = (pividx, a)
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
                    

# dummy methods if we don't want to trace
function new_tr_mat(nrows::Int,
                    tr::NoTracer)

    return NoTracerMatrix()
end

function add_row!(tr_mat::NoTracerMatrix,
                  sig::Sig,
                  row_ind::Int,
                  parent_ind::Int,
                  nred::Int)

    return
end

function resize_tracer_row_ops!(tr_mat::NoTracerMatrix,
                                row_ind::Int,
                                sz::Int)

    return
end

function store_row_op!(tr_mat::NoTracerMatrix,
                       row_ind::Int,
                       pividx::Int,
                       a::Coeff)

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
