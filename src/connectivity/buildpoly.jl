"""
    change_ringvar(F, new_S::Vector{Symbol})

Inject polynomial(s) into a new multivariate polynomial ring defined by `new_S`.
Uses the direct ring constructor `R(coeffs, exps)` to automatically handle
monomial ordering, coalescing, and zero-term removal.
"""
function change_ringvar(F::AbstractVector{<:Union{PolyRingElem, MPolyRingElem}}, new_S::Vector{Symbol})
    isempty(F) && return typeof(first(F))[]

    R = parent(first(F))
    old_S = symbols(R)
    newR, _ = polynomial_ring(base_ring(R), new_S)

    # Pre-calculate index mapping: old_index -> new_index (0 if removed)
    idx_map = [isnothing(idx) ? 0 : idx for idx in (findfirst(==(v), new_S) for v in old_S)]

    T = typeof(zero(base_ring(R)))
    n_vars = length(new_S)

    return map(F) do f
        coeffs = T[]
        exps = Vector{Int}[]

        if f isa MPolyRingElem
            for (e, c) in zip(exponent_vectors(f), coefficients(f))
                e_new = zeros(Int, n_vars)
                for (i, power) in enumerate(e)
                    if power > 0
                        @assert idx_map[i] != 0 "Occurrence of removed variable '$(old_S[i])' found!"
                        e_new[idx_map[i]] = power
                    end
                end
                push!(coeffs, c)
                push!(exps, e_new)
            end
        else # PolyRingElem
            @assert idx_map[1] != 0 "Occurrence of removed variable '$(old_S[1])' found!"
            for (deg_plus_1, c) in enumerate(coefficients(f))
                if !iszero(c)
                    e_new = zeros(Int, n_vars)
                    if deg_plus_1 > 1
                        e_new[idx_map[1]] = deg_plus_1 - 1
                    end
                    push!(coeffs, c)
                    push!(exps, e_new)
                end
            end
        end

        # Create the new polynomial associated to f
        return newR(coeffs, exps)
    end
end

# Dispatch for a single polynomial
change_ringvar(f::Union{PolyRingElem, MPolyRingElem}, new_S::Vector{Symbol}) = change_ringvar([f], new_S)[1]

 """
     change_ringvar(F)

Automatically detect all unique variables present in the polynomials and inject them into a ring containing only those variables.
"""
function change_ringvar(F::AbstractVector{<:MPolyRingElem})
    used_vars = unique(reduce(vcat, [Symbol.(vars(f)) for f in F]))
    return change_ringvar(F, used_vars)
end

change_ringvar(f::MPolyRingElem) = change_ringvar([f])[1]


"""
    MPolyBuild(F::Vector{Vector{RingElem}}, new_S::Vector{Symbol}, idx::Int)

Construct multivariate polynomials in a single variable.
`F` is a list of coefficient lists in degree-increasing order.
The polynomial will use the variable at index `idx` in `new_S`.
"""
function MPolyBuild(F::AbstractVector{<:AbstractVector{<:RingElement}}, new_S::Vector{Symbol}, idx::Int)
    isempty(F) && return []

    A = parent(first(first(F)))
    R, _ = polynomial_ring(A, new_S)

    T = typeof(zero(A))
    n_vars = length(new_S)

    return map(F) do f_coeffs
        coeffs = T[]
        exps = Vector{Int}[]

        # Construct coeff/exp lists
        for (deg_plus_1, c) in enumerate(f_coeffs)
            if !iszero(c)
                e_new = zeros(Int, n_vars)
                e_new[idx] = deg_plus_1 - 1
                push!(coeffs, c)
                push!(exps, e_new)
            end
        end

        # Create the polynomial associated to f_coeffs
        return R(coeffs, exps)
    end
end

# Dispatch for a single coefficient vector
MPolyBuild(f::AbstractVector{<:RingElement}, new_S::Vector{Symbol}, idx::Int) = MPolyBuild([f], new_S, idx)[1]

"""
    to_univariate(P::MPolyRingElem, var_idx::Int)

Converts a bivariate polynomial in `K[x, y]` into a univariate polynomial in `(K[x])[y]`.
"""
function to_univariate(P::MPolyRingElem, var_idx::Int)
    R = parent(P)
    K = base_ring(R)
    coeff_idx = 3 - var_idx # The variable to be pushed to the coefficients

    Rx, x = polynomial_ring(K, symbols(R)[coeff_idx])
    Rxy, y = polynomial_ring(Rx, symbols(R)[var_idx])

    res = zero(Rxy)
    for (c, exp) in zip(coefficients(P), exponent_vectors(P))
        # Build the coefficient in Rx, then attach to y in Rxy
        res += (c * x^exp[coeff_idx]) * y^exp[var_idx]
    end
    return res
end


"""
    to_bivariate(P::PolyRingElem, R_orig::MPolyRing)

Maps a univariate polynomial in `(K[x])[y]` back to the original bivariate ring `K[x, y]`.
"""
function to_bivariate(P::PolyRingElem, R_orig::MPolyRing)
    ctx = MPolyBuildCtx(R_orig)
    for (deg_y, c_y) in enumerate(coefficients(P))
        for (deg_x, c_x) in enumerate(coefficients(c_y))
            push_term!(ctx, c_x, [deg_x - 1, deg_y - 1])
        end
    end
    return finish(ctx)
end