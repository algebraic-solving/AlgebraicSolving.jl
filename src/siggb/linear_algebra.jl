function echelonize!(matrix::MacaulayMatrix,
                     tr::Tracer,
                     char::Val{Char},
                     shift::Val{Shift}) where {Char, Shift}

    arit_ops = 0

    col2hash = matrix.col2hash
    buffer = zeros(Cbuf, matrix.ncols)
    hash2col = Vector{MonIdx}(undef, matrix.ncols)
    rev_sigorder = Vector{Int}(undef, matrix.nrows)
    pivots = matrix.pivots

    tr_mat = TracerMatrix(Dict{Int, Int}(),
                          Coeff[],
                          Vector{Tuple{Int, Coeff}}[])

    @inbounds for i in 1:matrix.nrows
        rev_sigorder[matrix.sig_order[i]] = i
        row_ind = matrix.sig_order[i]
    end

    @inbounds for i in 1:matrix.ncols
        hash2col[col2hash[i]] = MonIdx(i)
    end

    n_trmat_rows = 0 
    @inbounds for i in 1:matrix.nrows
        row_ind = matrix.sig_order[i]

        row_cols = matrix.rows[row_ind]
        l_col_idx = hash2col[first(row_cols)]
        pivots[l_col_idx] == row_ind && continue

        # check if the row can be reduced
        does_red = false
        for (j, m_idx) in enumerate(row_cols)
            colidx = hash2col[m_idx]
            pividx = pivots[colidx]
            does_red = !iszero(pividx) && rev_sigorder[pividx] < i
            does_red && break
        end
        if !does_red
            pivots[l_col_idx] = row_ind
            continue
        end

        # indicate that we are tracing for row i
        n_trmat_rows += 1
        tr_mat.row_inds[row_ind] = n_trmat_rows

        # buffer the row
        row_coeffs = matrix.coeffs[row_ind]
        @inbounds for (k, j) in enumerate(row_cols)
            col_idx = hash2col[j]
            buffer[col_idx] = row_coeffs[k]
        end

        # allocate a row for the tracer matrix
        # at most we subtract (i-1) other rows
        row_ops = Vector{Tuple{Int, Coeff}}(undef, i - 1)
        
        # do the reduction
        n_row_subs = 0
        @inbounds for j in 1:matrix.ncols
            a = buffer[j] % Char
            iszero(a) && continue
            pividx = pivots[j]
            if iszero(pividx) || rev_sigorder[pividx] >= i
                continue
            end

            n_row_subs += 1
            row_ops[n_row_subs] = (pividx, a)

            # subtract a*rows[pivots[j]] from buffer
            pivmons = matrix.rows[pividx]
            pivcoeffs = matrix.coeffs[pividx]

            arit_ops_new = critical_loop!(buffer, j, a, hash2col, pivmons,
                                          pivcoeffs, shift)
            arit_ops += arit_ops_new
        end

        # finalize tracer row, add it to tracer matrix
        resize!(row_ops, n_row_subs)
        push!(tr_mat.col_inds_and_coeffs, row_ops)

        new_row_length = 0
        @inbounds for j in 1:matrix.ncols
            iszero(buffer[j]) && continue
            buffer[j] = buffer[j] % Char
            iszero(buffer[j]) && continue
            new_row_length += 1
        end

        # write out matrix row again
        j = 1
        inver = one(Coeff)
        new_row = Vector{MonIdx}(undef, new_row_length)
        new_coeffs = Vector{Coeff}(undef, new_row_length)
        @inbounds for k in 1:matrix.ncols
            iszero(buffer[k]) && continue
            if isone(j)
                pivots[k] = row_ind
                inver = inv(Coeff(buffer[k]), char)
            end
            new_row[j] = col2hash[k]
            new_coeffs[j] = isone(j) ? one(Coeff) : mul(inver, buffer[k], char)
            buffer[k] = zero(Cbuf)
            j += 1
        end
        # store that we normalized the row
        push!(tr_mat.diagonal, inver)

        # check if row lead reduced, TODO: dont know if this is reliable
        s = matrix.sigs[row_ind]
        m = monomial(s)
        @inbounds if isempty(new_row) || (matrix.rows[row_ind][1] != new_row[1] && any(!iszero, m.exps))
            matrix.toadd[matrix.toadd_length+1] = row_ind
            matrix.toadd_length += 1
        end

        matrix.rows[row_ind] = new_row
        matrix.coeffs[row_ind] = new_coeffs
    end
    if !iszero(arit_ops)
        @info "$(arit_ops) submul's"
    end

    # extend the tracer
    mat_basis_inds = Vector{Int}(undef, matrix.nrows)
    @inbounds for i in 1:matrix.nrows
        mat_basis_inds[i] = matrix.parent_inds[i]
    end
    push!(tr.basis_indices, mat_basis_inds)
    @inbounds push!(tr.sigs, matrix.sigs[1:matrix.nrows])
    push!(tr.matrices, tr_mat)

    return
end

# subtract mult
# TODO: for module tracking we won't be able to assume that mult = buffer[bufind]
@inline function critical_loop!(buffer::Vector{Cbuf},
                                bufind::Int,
                                mult::Cbuf,
                                hash2col::Vector{MonIdx},
                                pivmons::Vector{MonIdx},
                                pivcoeffs::Vector{Coeff},
                                shift::Val{Shift}) where Shift
    
    @inbounds buffer[bufind] = zero(Cbuf)
    l = length(pivmons)
    @turbo warn_check_args=false for k in 2:l
        c = pivcoeffs[k]
        m_idx = pivmons[k]
        colidx = hash2col[m_idx]
        buffer[colidx] = submul(buffer[colidx], mult, c, shift)
    end
    return l-1
end

# helper functions
# field arithmetic
function maxshift(::Val{Char}) where Char
    bufchar = Cbuf(Char)
    return bufchar << leading_zeros(bufchar)
end

# compute a representation of a - b*c mod char (char ~ Shift)
@inline function submul(a::Cbuf, b::Cbuf, c::Coeff, ::Val{Shift}) where Shift
    r0 = a - b*Cbuf(c)
    r1 = r0 + Shift
    r0 > a ? r1 : r0
end

@inline function inv(a::Coeff, ::Val{Char}) where Char
    return invmod(Cbuf(a), Cbuf(Char)) % Coeff
end

@inline function mul(a, b, ::Val{Char}) where Char 
    return Coeff((Cbuf(a) * Cbuf(b)) % Char)
end

@inline function add(a, b, ::Val{Char}) where Char
    c0 = a + b
    c1 = c0 - Coeff(Char)
    return max(c0, c1)
end

# for tracer

function new_tracer(::Val{N}) where N
    st_deg = 1
    basis_indices = Vector{Int}[]
    sigs = Vector{Sig{N}}[]
    mats = TracerMatrix[]
    return Tracer(st_deg, basis_indices, sigs, mats)
end

# for debug helping

function is_triangular(matrix::MacaulayMatrix)
    lms = [first(row) for row in matrix.rows[1:matrix.nrows] if !isempty(row)]
    return length(lms) == length(unique(lms))
end
