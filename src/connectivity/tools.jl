"""
    order_permut2d(L)

Return the indices that sort a 2D collection by its values.

Given a nested array `L` (e.g. a matrix or vector of vectors), this function
flattens all entries, sorts them by value, and returns the corresponding
2D indices.
"""
function order_permut2d(L)
    # Create a list of tuples with elements and their corresponding indices
    LL = [(L[i][j], (i, j)) for i in eachindex(L) for j in eachindex(L[i])]
    # Sort the enumerated list based on the values
    sorted_LL = sort(LL, by = x -> x[1])
    # Extract the sorted values and their corresponding indices
    sorted_ind = [pair[2] for pair in sorted_LL]

    return sorted_ind
end

"""
    trimat_rand(A, n; up=true, range=-100:100)

Generate a random triangular matrix over a given coefficient domain.
"""
function trimat_rand(A, n; up=true, range=-100:100)
    if up
        return [ i==j ? one(A) : (i<j ? A(rand(range)) : zero(A)) for i in 1:n, j in 1:n ]
    else
        return [ i==j ? one(A) : (i>j ? A(rand(range)) : zero(A)) for i in 1:n, j in 1:n ]
    end
end

"""
    int_coeffs(F::Vector{P}) where P <: Union{QQMPolyRingElem, QQPolyRingElem}
    int_coeffs(f::Union{QQMPolyRingElem, QQPolyRingElem})

Clear denominators of polynomial coefficients.

This function rescales each polynomial so that all coefficients become integers
(by multiplying with the least common multiple of denominators).
"""
function int_coeffs(F::Vector{P} where P <: Union{QQMPolyRingElem, QQPolyRingElem})
    CD = [ iszero(f) ? f : lcm(map(denominator, collect(coefficients(f)))) for f in F ]
    return (F .* CD)
end

int_coeffs(f::Union{QQMPolyRingElem, QQPolyRingElem}) = int_coeffs([f])[1]

macro iftime(v, ex)
    quote
        if $(esc(v))
            @time $(esc(ex))
        else
            $(esc(ex))
        end
    end
end


function param_use_lfs(I::Ideal, cfs_lfs::Vector{Vector{T}} where T<: RingElem, new_RS::Symbol)
    new_R = polynomial_ring(base_ring(I), new_RS)
    new_gens = [change_ringvar(f, new_RS) for f in gens(I)]
    I_new = Ideal(vcat(new_gens, [transpose(cfs_lf) * gens(new_R) for cfs_lf in cfs_lfs]))

    return rational_parametrization(I_new)
end