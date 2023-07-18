function echelonize!(matrix::MacaulayMatrix,
                     char::Val{Char},
                     shift::Val{Shift}) where {Char, Shift}


    col2hash = matrix.col2hash
    buffer = zeros(Cbuf, matrix.ncols)
    hash2col = Vector{MonIdx}(undef, matrix.ncols)
    rev_sigorder = Vector{Int}(undef, matrix.nrows)
    pivots = matrix.pivots

    @inbounds for i in 1:matrix.nrows
        rev_sigorder[matrix.sig_order[i]] = i
        row_ind = matrix.sig_order[i]
    end

    @inbounds for i in 1:matrix.ncols
        hash2col[col2hash[i]] = MonIdx(i)
    end

    @inbounds for i in 2:matrix.nrows
        row_ind = matrix.sig_order[i]

        row_cols = matrix.rows[row_ind]

        # check if the row can be reduced
        does_red = false
        for (j, m_idx) in enumerate(row_cols)
            colidx = hash2col[m_idx]
            pividx = pivots[colidx]
            does_red = !iszero(pividx) && rev_sigorder[pividx] < i
            does_red && break
        end
        !does_red && continue

        # buffer the row
        row_coeffs = matrix.coeffs[row_ind]
        @inbounds for (k, j) in enumerate(row_cols)
            col_idx = hash2col[j]
            buffer[col_idx] = row_coeffs[k]
        end

        # do the reduction
        @inbounds for j in 1:matrix.ncols
            a = buffer[j] % Char
            iszero(a) && continue
            pividx = pivots[j]
            if iszero(pividx) || rev_sigorder[pividx] >= i
                continue
            end

            # subtract m*rows[pivots[j]] from buffer
            pivcoeffs = matrix.coeffs[pividx]
            b = inv(pivcoeffs[1], char)
            m = mul(a, b, char)

            buffer[j] = zero(Cbuf)
            @inbounds for (k, m_idx) in enumerate(matrix.rows[pividx])
                isone(k) && continue
                c = pivcoeffs[k]
                colidx = hash2col[m_idx]
                buffer[colidx] = submul(buffer[colidx], m, c, shift)
            end
        end

        new_row_length = 0
        @inbounds for j in 1:matrix.ncols
            buffer[j] = buffer[j] % Char
            iszero(buffer[j]) && continue
            new_row_length += 1
        end

        # write out matrix row again
        j = 1
        new_row = Vector{MonIdx}(undef, new_row_length)
        new_coeffs = Vector{Coeff}(undef, new_row_length)
        @inbounds for k in 1:matrix.ncols
            iszero(buffer[k]) && continue
            new_row[j] = col2hash[k]
            new_coeffs[j] = buffer[k]
            if isone(j)
                pivots[k] = row_ind
            end
            buffer[k] = zero(Cbuf)
            j += 1
        end

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
end


# helper functions
# field arithmetic
function maxshift(::Val{Char}) where Char
    bufchar = Cbuf(Char)
    return bufchar << leading_zeros(bufchar)
end

# compute a representation of a - b*c mod char (char ~ Shift)
@inline function submul(a::Cbuf, b::Coeff, c::Coeff, ::Val{Shift}) where Shift
    r0 = a - Cbuf(b)*Cbuf(c)
    r1 = r0 + Shift
    r0 > a ? r1 : r0
end

@inline function inv(a::Coeff, ::Val{Char}) where Char
    return invmod(Cbuf(a), Cbuf(Char)) % Coeff
end

@inline function mul(a, b, ::Val{Char}) where Char 
    return Coeff((Cbuf(a) * Cbuf(b)) % Char)
end
