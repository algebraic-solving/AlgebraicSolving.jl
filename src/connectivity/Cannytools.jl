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

    ## Is it correct AND useful?
    sort!(J)
    ##
    Jphi = [ j for j in J if j <= nphi ]
    Jproj = setdiff(J, Jphi)

    # Construct the truncated Jacobian matrix
    psi = vcat(V.gens, phi[Jphi])
    JW = matrix(R, QQMPolyRingElem[ derivative(f, k) for f in psi, k in setdiff(1:n, Jproj)])
    # Compute the minors
    sizeminors = c + length(Jphi) + min(dimproj, length(J)-1) - (length(J)-1)
    minors = _compute_minors(sizeminors, JW, R)

    if only_mins
        return minors
    else
        return vcat(V.gens, minors)
    end
end

function _compute_minors(p, A, R)
    #Computes the p-minors of a matrix A
    n, m = size(A)
    rowsmins = collect(_combinations(1:n, p, 1, Vector{Int}([])))
    colsmins = collect(_combinations(1:m, p, 1, Vector{Int}([])))
    mins = Vector{eltype(A)}(undef, length(rowsmins) * length(colsmins))
    k = 1
    # for performance tweaks, check if there are non-linear polynomials or if all are linear or constant
    degmax = 0
    for a in A
        degmax = min(2, max(degmax, total_degree(a)))
        if degmax > 1
            break
        end
    end
    # TODO: use algos from FLINT (in particular berkowitz)
    # Naive function performs better for high degree or small matrices
    if (degmax == 0 && p <= 2) || (degmax == 1 && p <= 6) || (degmax == 2)
        detfct = s->detmpoly(s, R)
    else # else use fraction-free LU from AbstractAlgebra.jl
        detfct = det
    end
    for rowsmin in rowsmins
        for colsmin in colsmins
            mins[k] = detfct(A[rowsmin, colsmin])
            k += 1
        end
    end

    return mins
end

function _combinations(a, n, start, chosen)
    if length(chosen) == n
        return [chosen]
    elseif start > length(a)
        return Vector{Int}([])
    else
        # Include the current element and recurse
        include_current = _combinations(a, n, start + 1, [chosen; a[start]])
        # Exclude the current element and recurse
        exclude_current = _combinations(a, n, start + 1, chosen)
        return vcat(include_current, exclude_current)
    end
end

function detmpoly(A, R)
    # Get the size of the matrix
    n = size(A, 1)
    if n != size(A, 2)
        throw(ArgumentError("Matrix must be square"))
    end

    if n == 1
        return A[1, 1]
    end

    if n == 2
        return A[1, 1] * A[2, 2] - A[1,2] * A[2, 1]
    end

    # Initialize the determinant polynomial
    detA = zero(R)

    # Compute the determinant polynomial
    for j = 1:n
        submatrix = A[2:end, [i for i = 1:n if i != j]]
        detA += (-1)^(1+j)*A[1, j] * detmpoly(submatrix, R)
    end

    return detA
end

function MidRationalPoints(S::Vector{Vector{T}}, Q::Vector{T} = T[]) where {T <: QQFieldElem}
    # * S is a list of [ [l_1,r_1], ..., [l_n, r_n] ]
    # such that the [l_i, r_i] are rational and disjoint open intervals.
    # * Q is a list of rationals that intersects no [l_i,r_i]
    #
    # It orders the [l_i,r_i], and compute a list ratioP such that
    # strictly between each of these intervals there is:
    # - either at least one element of Q
    # - or the simplest rational number
    isempty(S) && return Q

    S1, Q1 = sort(S, lt=(x, y) -> x[2] <= y[1]), sort(Q)
    ratioP = T[]
    qidx = 1
    qlen = length(Q1)

    # Handle left gap before first interval
    while qidx <= qlen && Q1[qidx] < S1[1][1]
        push!(ratioP, Q1[qidx])
        qidx += 1
    end

    # Loop through gaps between sorted disjoint intervals
    for i in 1:(length(S1) - 1)
        ri, li1 = S1[i][2], S1[i+1][1]
        @assert ri < li1 "Intervals are not disjoint."
        inserted = false
        while qidx <= qlen && Q1[qidx] < li1
            @assert(Q1[qidx] > ri, "A query point is singular")
            push!(ratioP, Q1[qidx])
            inserted = true
            qidx += 1
        end
        if !inserted
            eps = (li1 - ri)//1000 # for open interval
            # We choose the simplest in absolute value
            if -ri > li1 # this means ri is negative and the largest in absolute value
                push!(ratioP, -simplest_between(-ri - eps, -li1 + eps))
            else
                push!(ratioP, simplest_between(ri + eps, li1 - eps))
            end
        end
    end

    # Append remaining right-side Q points
    while qidx <= qlen
        push!(ratioP, Q1[qidx])
        qidx += 1
    end

    return ratioP
end

