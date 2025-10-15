# Compute the set of points x where, if psi in the map whose elts are the ones in [phi_1,..,phi_p, x_{p+1},...x_n] indexed in J,
# then psi(T_x(V)) has dimension < dimproj
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
    JW = transpose([ derivative(f, k) for k in setdiff(1:n, Jproj), f in psi])
    # Compute the minors
    sizeminors = c + length(Jphi) + min(dimproj, length(J)-1) - (length(J)-1)
    minors = compute_minors(sizeminors, JW, R)

    if only_mins
        return minors
    else
        return vcat(V.gens, minors)
    end
end

function compute_minors(p, A, R)
    #Computes the p-minors of a matrix A
    n, m = size(A)
    rowsmins = collect(combinations(1:n, p))
    colsmins = collect(combinations(1:m, p))
    mins = Vector{eltype(A)}(undef, length(rowsmins) * length(colsmins))
    k = 1
    for rowsmin in rowsmins
        for colsmin in colsmins
            mins[k] = detmpoly(A[rowsmin, colsmin], R)
            k += 1
        end
    end

    return mins
end

function combinations(a, n)
    # Helper function to recursively generate combinations
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
    return _combinations(a, n, 1, Vector{Int}([]))
end

function detmpoly(A::Matrix{T} where T<:MPolyRingElem, R)
    # Get the size of the matrix
    n = size(A, 1)
    if n != size(A, 2)
        throw(ArgumentError("Matrix must be square"))
    end

    if n == 1
        return A[1, 1]
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

