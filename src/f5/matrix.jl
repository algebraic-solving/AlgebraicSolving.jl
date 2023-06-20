# TODO: update this function
function initialize_matrix(::Type{C}) where {C}
    rows = Vector{Vector{ColIdx}}(undef, 0)
    col2hash = Vector{MonIdx}(undef, 0)
    coeffs = Vector{Vector{C}}(undef, 0)

    size = 0
    npivots = 0
    nrows = 0
    ncols = 0

    MacaulayMatrix(rows, col2hash, coeffs, size, npivots,
                   nrows, ncols)
end

# Refresh and initialize matrix for `npairs` elements
function reinitialize_matrix!(matrix::MacaulayMatrix, npairs::Int)
    matrix.size = 2 * npairs
    resize!(matrix.rows, matrix.size)
    resize!(matrix.sig_order, matrix.size)
    resize!(matrix.coeffs, matrix.size)
    resize!(matrix.pivots, 2 * npairs)
    matrix.pivot_size = 2 * npairs
    return matrix
end
