@doc Markdown.doc"""
    roadmap(I::Ideal{T} where T <: MPolyRingElem, <keyword arguments>)

Given a **radical** ideal `I` with solution set X, that is smooth and
whose real trace XR is bounded, return a roadmap of XR

The output is given as a Roadmap structure, encoding the recursive structure
of roadmaps. It is encoded as a chained list, whose root is containing the equations defining X
and each node representing a curve component, that is defined by additional polar equation and base point.
Moreover it is linked to fibers, that share the same base point.

# Arguments
- `I::Ideal{T} where T <: QQMPolyRingElem`: input generators.
- `C::Vector{Vector{QQFieldElem}}=Vector{QQFieldElem}[]`: query points with rational coefficients
- `info_level::Int=0`: verbosity level
- `checks::Bool=false`: whether perform checks (dimension, regularity, etc.)
- 'generic_change=false": whether it performs a prior random linear change of variables (TODO)
)

# Examples
```jldoctest
julia> using AlgebraicSolving

julia> R,(x1,x2,x3,x4) = polynomial_ring(QQ, ["x1","x2","x3","x4"])
(Multivariate polynomial ring in 4 variables over QQ, QQMPolyRingElem[x1, x2, x3, x4])

julia> I = Ideal([(x1^2+x2^2+x3^2+x4^2+9-1)^2-4*9*(x1^2+x2^2+x3^2) + 1])
QQMPolyRingElem[x1^4 + 2*x1^2*x2^2 + 2*x1^2*x3^2 + 2*x1^2*x4^2 - 20*x1^2 + x2^4 + 2*x2^2*x3^2 + 2*x2^2*x4^2 - 20*x2^2 + x3^4 + 2*x3^2*x4^2 - 20*x3^2 + x4^4 + 16*x4^2 + 65]

julia> RM = roadmap(I, checks=true)
Vector{QQFieldElem}[[], [-3], [-3, -2], [-2], [-2, -1], [-2, 0], [-2, 1], [3], [3, 2]]

julia> nb_nodes(RM)
9

julia> all_eqs(RM)
9-element Vector{Ideal{QQMPolyRingElem}}:
 QQMPolyRingElem[x1^4 + 2*x1^2*x2^2 + 2*x1^2*x3^2 + 2*x1^2*x4^2 - 20*x1^2 + x2^4 + 2*x2^2*x3^2 + 2*x2^2*x4^2 - 20*x2^2 + x3^4 + 2*x3^2*x4^2 - 20*x3^2 + x4^4 + 16*x4^2 + 65, -4*x1^2*x3 - 4*x2^2*x3 - 4*x3^3 - 4*x3*x4^2 + 40*x3, -4*x1^2*x4 - 4*x2^2*x4 - 4*x3^2*x4 - 4*x4^3 - 32*x4]
 QQMPolyRingElem[x1^4 + 2*x1^2*x2^2 + 2*x1^2*x3^2 + 2*x1^2*x4^2 - 20*x1^2 + x2^4 + 2*x2^2*x3^2 + 2*x2^2*x4^2 - 20*x2^2 + x3^4 + 2*x3^2*x4^2 - 20*x3^2 + x4^4 + 16*x4^2 + 65, -4*x1^2*x4 - 4*x2^2*x4 - 4*x3^2*x4 - 4*x4^3 - 32*x4, x1 + 3]
 QQMPolyRingElem[x1^4 + 2*x1^2*x2^2 + 2*x1^2*x3^2 + 2*x1^2*x4^2 - 20*x1^2 + x2^4 + 2*x2^2*x3^2 + 2*x2^2*x4^2 - 20*x2^2 + x3^4 + 2*x3^2*x4^2 - 20*x3^2 + x4^4 + 16*x4^2 + 65, x1 + 3, x2 + 2]
 QQMPolyRingElem[x1^4 + 2*x1^2*x2^2 + 2*x1^2*x3^2 + 2*x1^2*x4^2 - 20*x1^2 + x2^4 + 2*x2^2*x3^2 + 2*x2^2*x4^2 - 20*x2^2 + x3^4 + 2*x3^2*x4^2 - 20*x3^2 + x4^4 + 16*x4^2 + 65, -4*x1^2*x4 - 4*x2^2*x4 - 4*x3^2*x4 - 4*x4^3 - 32*x4, x1 + 2]
 QQMPolyRingElem[x1^4 + 2*x1^2*x2^2 + 2*x1^2*x3^2 + 2*x1^2*x4^2 - 20*x1^2 + x2^4 + 2*x2^2*x3^2 + 2*x2^2*x4^2 - 20*x2^2 + x3^4 + 2*x3^2*x4^2 - 20*x3^2 + x4^4 + 16*x4^2 + 65, x1 + 2, x2 + 1]
 QQMPolyRingElem[x1^4 + 2*x1^2*x2^2 + 2*x1^2*x3^2 + 2*x1^2*x4^2 - 20*x1^2 + x2^4 + 2*x2^2*x3^2 + 2*x2^2*x4^2 - 20*x2^2 + x3^4 + 2*x3^2*x4^2 - 20*x3^2 + x4^4 + 16*x4^2 + 65, x1 + 2, x2]
 QQMPolyRingElem[x1^4 + 2*x1^2*x2^2 + 2*x1^2*x3^2 + 2*x1^2*x4^2 - 20*x1^2 + x2^4 + 2*x2^2*x3^2 + 2*x2^2*x4^2 - 20*x2^2 + x3^4 + 2*x3^2*x4^2 - 20*x3^2 + x4^4 + 16*x4^2 + 65, x1 + 2, x2 - 1]
 QQMPolyRingElem[x1^4 + 2*x1^2*x2^2 + 2*x1^2*x3^2 + 2*x1^2*x4^2 - 20*x1^2 + x2^4 + 2*x2^2*x3^2 + 2*x2^2*x4^2 - 20*x2^2 + x3^4 + 2*x3^2*x4^2 - 20*x3^2 + x4^4 + 16*x4^2 + 65, -4*x1^2*x4 - 4*x2^2*x4 - 4*x3^2*x4 - 4*x4^3 - 32*x4, x1 - 3]
 QQMPolyRingElem[x1^4 + 2*x1^2*x2^2 + 2*x1^2*x3^2 + 2*x1^2*x4^2 - 20*x1^2 + x2^4 + 2*x2^2*x3^2 + 2*x2^2*x4^2 - 20*x2^2 + x3^4 + 2*x3^2*x4^2 - 20*x3^2 + x4^4 + 16*x4^2 + 65, x1 - 3, x2 - 2]
```
"""
function roadmap(
    I::Ideal{P};                                                            # input ideal
    C::Vector{Vector{Vector{QQFieldElem}}}=Vector{Vector{QQFieldElem}}[],   # query points: interval with rational coefficients
    info_level::Int=0,                                                      # verbosity level
    checks::Bool=false                                                      # perform checks (dimension, regularity, etc.)
) where {P <: QQMPolyRingElem}
    # L_chosen holds the sequence of generic linear forms shared across ALL fibers.
    L_chosen = P[]
    empty_base_pt = Tuple{P, QQFieldElem}[]

    return _roadmap_rec(I, empty_base_pt, C, L_chosen, info_level, checks)
end

@doc Markdown.doc"""
    roadmap(I::Ideal{T} where T <: MPolyRingElem, I::Ideal{P}, C::Ideal{P}; info_level::Int=0, checks::Bool=false)
```
"""
function roadmap(
    I::Ideal{P},                # input ideal
    C::Ideal{P};                # ideal defining query points
    info_level::Int=0,          # verbosity level
    checks::Bool=false          # perform checks (dimension, regularity, etc.)
) where {P <: QQMPolyRingElem}
    @assert(parent(I)==parent(C), "Equations for variety and query points must live the same ring")
    CI = real_solutions(C, info_level=max(info_level-1,0), nr_thrds=Threads.nthreads())
    return roadmap(I, CI, info_level, checks)
end

function _roadmap_rec(
    I::Ideal{P},     # input ideal
    base_pt::BasePoint{P},                      # single base point with rational coefficients
    C::Vector{Vector{Vector{QQFieldElem}}},     # query points with rational coefficients
    L_chosen::Vector{P},
    info_level::Int,                            # verbosity level
    checks::Bool                                # perform checks (dimension, regularity, etc.)
) where {P <: QQMPolyRingElem}

    n = nvars(parent(I))
    if n <= 2
        return Roadmap(I, RMnode(base_pt, I.gens, RMnode[]))
    end

    # Some preprocessng
    if isnothing(I.dim)
        lucky_prime = _generate_lucky_primes(I.gens, one(ZZ)<<30, one(ZZ)<<31-1, 1) |> first
        INEW = Ideal(change_base_ring.(Ref(GF(lucky_prime)), I.gens))
        I.dim = dimension(INEW)
    end

    e = length(base_pt)

    # Base case (dim(I) <= 1)
    if I.dim <= 1
        return Roadmap(I, RMnode(base_pt, I.gens, RMnode[]))
    end
    # Terminal case (dim(F) <= 1)
    if I.dim - e <= 1
        return RMnode(base_pt, P[], RMnode{P}[])
    end

    # If we reach a new depth level, then generate new linear forms
    if length(L_chosen) < e + 2
        info_level > 0 && println("Generating generic linear forms for depth $e...")
        #append_generic_linear_forms!(L_chosen, I, base_pt, info_level, checks)
        if e == 0
            push!(L_chosen, gens(parent(I))[1])
        end
        push!(L_chosen, gens(parent(I))[e+2])
    end
    L_e1 = L_chosen[1 : e+1]
    L_e2 = L_chosen[1 : e+2]
    L_next = L_chosen[e+1] # The linear form we step along in this iteration

    ## K(L_e1,Fq) ##
    info_level > 0 && println("Compute L$(e+1)-critical points: K1")
    K1Fq_pol = computepolar(1:e+1, I, phi=L_e1) |> Ideal # TODO: reuse for different fiber but same e
    K1Fq_fiber = vcat(K1Fq_pol.gens, [l - q for (l, q) in base_pt]) |> Ideal
    K1Fq = real_solutions(K1Fq_fiber, info_level=max(info_level-1,0), nr_thrds=Threads.nthreads(), interval=true)

    ## W = K(L_e2, Fq) ##
    info_level>0 && println("Compute L_e2-polar variety: W")
    K2Fqmins = computepolar(1:e+2, I, phi=L_e2, only_mins=true)
    K2Fq = vcat(I.gens, K2Fqmins) |> Ideal
    K2Fq.dim = e + 1
    RM = RMnode(base_pt, K2Fqmins, RMnode{P}[])

    ## Points with vertical tangents in W ##
    info_level>0 && println("Compute W-critical points with <<vertical>> tangent: K1W")
    K1WmFq_pol = computepolar(1:e+2, K2Fq, phi=L_e2, dimproj=e) |> Ideal
    K1WmFq_fiber = Ideal(vcat(K1WmFq_pol.gens, [l - q for (l, q) in base_pt]))
    K1WmFq = real_solutions(K1WmFq_fiber, info_level=max(info_level-1,0), nr_thrds=Threads.nthreads(), interval=true)

    ## New base and query points ##
    # Cq = isempty(q) ? C : [ c for c in C if q[e] in c[e]]
    Cq = isempty(base_pt) ? C : [c for c in C if base_pt[e][2] in _linear_interval_eval(L_chosen[e], c)]

    K1W = vcat(K1Fq, K1WmFq)

    # Heuristic to be proven (Reeb's th)
    #K1W = K1W[2:end-1]
    ##########

    # Evaluate L_{e+1} on the n-dimensional points to project them down to the line
    K1W_proj = [_linear_interval_eval(L_next, pt) for pt in K1W]
    Cq_proj  = unique([_linear_interval_eval(L_next, c) for c in Cq])

    K1WRat = _mid_rational_points_inter(K1W_proj, Cq_proj)

    # Construct the generalized base points for the children
    newQ = [vcat(base_pt, [(L_next, r)]) for r in K1WRat]

    # --- 6. Recurse ---
    if !isempty(newQ)
        for newq in newQ
            RMFq = _roadmap_rec(I, newq, Cq, L_chosen, info_level, checks)
            push!(RM.children, RMFq)
        end
    end

    return e == 0 ? Roadmap(I, RM) : RM
end


function _mid_rational_points_inter(S::Vector{Vector{T}}, Q::Vector{Vector{T}} = Vector{T}[]) where {T <: QQFieldElem}
    # * S is a list of [ [l_1,r_1], ..., [l_n, r_n] ]
    # such that the [l_i, r_i] are rational and disjoint open intervals.
    # * Same assumptions on Q
    # * Intervals in S and Q do not intersect as well
    #
    # It orders the [l_i,r_i], and compute a list ratioP
    # intersecting all intervals of Q and such that
    # *strictly* between each of the [l_i,r_i] there is:
    # - either at least one element inside an interval of Q
    # - or the simplest rational number
    isempty(S) && return Q

    S1, Q1 = sort(S, lt=(x, y) -> x[2] <= y[1]), sort(Q, lt=(x, y) -> x[2] <= y[1])
    ratioP = T[]
    qidx = 1
    qlen = length(Q1)

    # Handle left gap before first interval
    while qidx <= qlen && Q1[qidx][2] < S1[1][1]
        ql, qr = Q1[qidx]
        push!(ratioP, _open_simplest_between(ql, qr, abs(qr - ql)//1000))
        qidx += 1
    end

    # Loop through gaps between sorted disjoint intervals
    for i in 1:(length(S1) - 1)
        ri, li1 = S1[i][2], S1[i+1][1]
        @assert ri < li1 "Intervals are not disjoint."
        inserted = false
        while qidx <= qlen && Q1[qidx][2] < li1
            ql, qr = Q1[qidx]
            @assert(ql > ri, "A query point might be singular")
            push!(ratioP, _open_simplest_between(ql, qr, abs(qr - ql)//1000))
            inserted = true
            qidx += 1
        end
        @assert qidx > qlen || Q1[qidx][1] > S1[i+1][2] "A query point might be singular"
        # If there's already rational between no need to add new
        !inserted && push!(ratioP, _open_simplest_between(ri, li1, abs(li1 - ri)//1000))
    end

    # Append remaining right-side Q points
    while qidx <= qlen
        ql, qr = Q1[qidx]
        push!(ratioP, _open_simplest_between(ql, qr, abs(qr - ql)//1000))
        qidx += 1
    end

    return ratioP
end

function _open_simplest_between(a::QQFieldElem, b::QQFieldElem, eps::QQFieldElem)
    # We choose the simplest in absolute value
    if -a > b # this means a is negative and the largest in absolute value
        return -simplest_between(-a - eps, -b + eps)
    else
        return  simplest_between( a + eps,  b - eps)
    end
end



function _linear_interval_eval(L::QQMPolyRingElem,
                       I::Vector{Vector{QQFieldElem}})
    R = parent(L)
    x = gens(R)

    lo = hi = zero(QQ)

    for (xi, (a, b)) in zip(x, I)
        c = coeff(L, xi)
        if c >= 0
            lo += c*a
            hi += c*b
        else
            lo += c*b
            hi += c*a
        end
    end

    return [lo, hi]
end