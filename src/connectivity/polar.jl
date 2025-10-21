@doc Markdown.doc"""
    computepolar(
        J::Union{Vector{Int},UnitRange{Int}},
        V::Ideal{P};
        phi::Vector{P}=P[],
        dimproj::Int=length(J)-1,
        only_mins::Bool=false
    ) where {P <: MPolyRingElem}

Compute the **polar variety** associated with the map whose components are
`[phi_1, …, phi_p, x_{p+1}, …, x_n]` indexed by `J`, for an algebraic variety defined
by the ideal `V`.

More precisely, this function computes the set of points `x` in `V` such that,
if `psi` denotes the above map, then the image of the tangent space `Tₓ(V)` under
the differential `psi` has dimension strictly less than `dimproj`.

This is a key geometric construction in the computation of polar varieties,
used in the roadmap algorithm and critical point methods.

# Arguments
- `J::Union{Vector{Int},UnitRange{Int}}`:
  Indices of the selected coordinate functions
- `V::Ideal{P} where P <: MPolyRingElem`:
  Input ideal defining the variety V on which critical loci is computed
- `phi::Vector{P}=P[]`:
  Polynomial map possibly completed by projection coordinates to have a map of total length `n`.
- `dimproj::Int=length(J)-1`:
  Expected maximal dimension of the image of the tangent space.
  Typically equals the number of projection coordinates minus one.
- `only_mins::Bool=false`:
  If `true`, only the computed minors of the Jacobian are returned;
  otherwise, the output includes both the generators of `V` and those minors.

# Returns
- If `only_mins=false` (default):
  A `Vector{P}` containing the union of the generators of `V` and the computed minors.
  This defines the ideal of the polar variety.
- If `only_mins=true`:
  A `Vector{P}` containing only the minors (without the original equations of `V`).

# Example
```jldoctest
julia> using AlgebraicSolving

julia> R,(x1,x2,x3,x4) = polynomial_ring(QQ, ["x1","x2","x3","x4"])
(Multivariate polynomial ring in 4 variables over QQ, QQMPolyRingElem[x1, x2, x3, x4])

julia> I = Ideal([(x1^2+x2^2+x3^2+x4^2+9-1)^2-4*9*(x1^2+x2^2+x3^2)])
QQMPolyRingElem[x1^4 + 2*x1^2*x2^2 + 2*x1^2*x3^2 + 2*x1^2*x4^2 - 20*x1^2 + x2^4 + 2*x2^2*x3^2 + 2*x2^2*x4^2 - 20*x2^2 + x3^4 + 2*x3^2*x4^2 - 20*x3^2 + x4^4 + 16*x4^2 + 64]

julia> computepolar(1:3, I, dimproj=1, phi=[x1^2+x2^2+x3^2+x4^2])
5-element Vector{QQMPolyRingElem}:
 x1^4 + 2*x1^2*x2^2 + 2*x1^2*x3^2 + 2*x1^2*x4^2 - 20*x1^2 + x2^4 + 2*x2^2*x3^2 + 2*x2^2*x4^2 - 20*x2^2 + x3^4 + 2*x3^2*x4^2 - 20*x3^2 + x4^4 + 16*x4^2 + 64
 4*x1^3 + 4*x1*x2^2 + 4*x1*x3^2 + 4*x1*x4^2 - 40*x1
 4*x1^2*x4 + 4*x2^2*x4 + 4*x3^2*x4 + 4*x4^3 + 32*x4
 2*x1
 2*x4
```
"""
function computepolar(
        J::Union{Vector{Int},UnitRange{Int}},   # coordinate images of [phi,proj] (as above)
        V::Ideal{P};                            # input ideal
        phi::Vector{P} = P[],                   # polynomial map in consideration (completed by sufficiently many projections)
        dimproj = length(J)-1,                  # maximum dimension of tangent space of phi
        only_mins = false                       # return only minors without eqns of V
    ) where (P <: MPolyRingElem)

    R = parent(V)
    n = nvars(R)
    nphi = length(phi)
    @assert all([1<=j<=n for j in J])

    isnothing(V.dim) && dimension(V)
    c = n - V.dim

    ##
    Jphi = [ j for j in J if j <= nphi ]
    Jproj = setdiff(J, Jphi)

    # Construct the truncated Jacobian matrix
    psi = vcat(V.gens, phi[Jphi])
    JW = matrix(R, QQMPolyRingElem[ derivative(f, k) for f in psi, k in setdiff(1:n, Jproj)])
    # Compute the minors
    sizeminors = c + length(Jphi) + min(dimproj, length(J)-1) - (length(J)-1)
    minors = _compute_minors(sizeminors, JW)

    if only_mins
        return minors
    else
        return vcat(V.gens, minors)
    end
end

function _compute_minors(p, A)
    # Computes the p-minors of a matrix A
    rowsmins = _combinations(1:nrows(A), p)
    colsmins = _combinations(1:ncols(A), p)
    # We use charpoly for a division-free determinant method
    return [ coeff(charpoly(A[rows, cols]), 0) for rows in rowsmins for cols in colsmins ]
end

function _combinations(v::UnitRange{Int}, k::Int) 
    # Compute the k-subsets of v
    n = length(v)
    ans = Vector{Int}[]
    k > n && return ans
    _combinations_dfs!(ans, Vector{Int}(undef, k), v, n, k)
    return ans
end

function _combinations_dfs!(ans::Vector{Vector{Int}}, comb::Vector{Int}, v::UnitRange{Int}, n::Int, k::Int)
    k < 1 && (pushfirst!(ans, comb[:]); return)
    for m in n:-1:k
        comb[k] = v[m]
        _combinations_dfs!(ans, comb, v, m - 1, k - 1)
    end
end