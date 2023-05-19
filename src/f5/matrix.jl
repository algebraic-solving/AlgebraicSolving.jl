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
    resize!(matrix.rows, npairs * 2)
    matrix.size = 2 * npairs
    matrix.ncols = 0
    matrix
end
