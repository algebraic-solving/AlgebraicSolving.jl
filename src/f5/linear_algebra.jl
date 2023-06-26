function buffer!(matrix::MacaulayMatrix{C},
                 buffer::Vector{C},
                 row_ind::Int) where C

    @inbounds row_cols = matrix.rows[row_ind]
    @inbounds row_coeffs = matrix.coeffs[row_ind]

    j = 1
    @inbounds for i in row_cols
        col_idx = matrix.hash2col[i]
        buffer[col_idx] = row_coeffs[j]
    end
    return
end
