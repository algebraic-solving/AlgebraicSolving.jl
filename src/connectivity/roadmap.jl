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
        append_generic_linear_forms!(L_chosen, I, base_pt)
    end
    L_e1 = L_chosen[1 : e+1]
    L_e2 = L_chosen[1 : e+2]
    L_next = L_chosen[e+1] # The linear form we step along in this iteration

    # TODO: add a variable for L_next and compute param w.r.t. it. Isolate only eliminating poly.
    # this directly garantees disjoint interval and correct enclosing
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

    if !isempty(newQ)
        for newq in newQ
            RMFq = _roadmap_rec(I, newq, Cq, L_chosen, info_level, checks)
            push!(RM.children, RMFq)
        end
    end

    return e == 0 ? Roadmap(I, RM) : RM
end


# -------------------------------------------------------------------------
# Genericity Tests & Generation (Completely Extracted)
# -------------------------------------------------------------------------

"""
    append_generic_linear_forms!(L_chosen, I, base_pt, info_level, checks)

Deterministically searches for and appends one linear form to `L_chosen` (two if empty)
so that all genericity conditions pass for the current fiber.
"""
function append_generic_linear_forms!(
    L_chosen::Vector{P},
    I::Ideal{P},
    base_pt::BasePoint{P},
) where {P <: QQMPolyRingElem}

    R = parent(I.gens[1])
    n = nvars(R)
    v = gens(R)

    # Simple deterministic search for coefficients
    idx = 1
    while true
        # Generate nb candidates deterministically based on idx
        nb = isempty(L_chosen) ? 2 : 1

        # TODO: Coefficient generation
        cfs = [[convert(Int64, hash((idx, j, i + length(L_chosen))) % 7) - 3 for j in 1:n] for i in 1:nb]
        @show cfs
        L_next = [ sum(cfs[i][j] * v[j] for j in 1:n) for i in 1:nb ]
        ####

        L_test = vcat(L_chosen, L_next)
        if check_genericity(I, base_pt, L_test)
            append!(L_chosen, L_next)
            return
        end
        idx += 1
    end
end

"""
    check_genericity(I, base_pt, L_test, checks)

Validates all generic roadmap assumptions
Returns true if `L_test` satisfies the geometric requirements.
"""
function check_genericity(
    I::Ideal{P},
    base_pt::BasePoint{P},
    L_test::Vector{P},
) where {P <: QQMPolyRingElem}
    lucky_prime = first(_generate_lucky_primes(I.gens, one(ZZ)<<30, (one(ZZ)<<31)-1, 1))
    Kp = GF(lucky_prime)
    Ip = Ideal(change_base_ring.(Ref(Kp), I.gens))
    Lp = change_base_ring.(Ref(Kp), L_test)
    fib_p = [ change_base_ring(Kp, l) - Kp(q) for (l, q) in base_pt]

    e = length(base_pt)
    Lp_e  = Lp[1 : e]
    Lp_e1 = Lp[1 : e+1]
    Lp_e2 = Lp[1 : e+2]
    Ip.dim = I.dim

    function dimfib(J)
        return dimension(Ideal(vcat(J, fib_p)))
    end

    # Test 1: Fiber dimension
    if dimfib(Ip.gens) != I.dim - e
        return false
    end

    # Test 2: Fiber must be smooth
    polar_e = computepolar(1:e, Ip, phi=Lp_e)
    if dimfib(polar_e) != -1
        return false
    end

    # Test 3: K1Fq must be finite
    K1Fq_pol = computepolar(1:e+1, Ip, phi=Lp_e1)
    if dimfib(K1Fq_pol) > 0
        return false
    end

    # Test 4: K2Fq must have dimension 1
    K2Fq_gens = computepolar(1:e+2, Ip, phi=Lp_e2)
    if dimfib(K2Fq_gens) != 1
        return false
    end

    # Test 5: Points with vertical tangents in W must be finite
    K1Wm = computepolar(1:e+2, Ideal(K2Fq_gens), phi=Lp_e2, dimproj=e)
    if dimfib(K1Wm) > 0
        return false
    end

    return true
end

# Generate N random primes between low and up
# that do not divide any numerator/denominator
# of any coefficient in polynomials from LP
function _generate_lucky_primes(
    LF::Vector{<:MPolyRingElem},
    low::ZZRingElem,
    up::ZZRingElem,
    N::Int64
    )
    # Using a Set avoids resizing and `unique!` shifting overhead
    CF_set = Set{ZZRingElem}()
    for f in LF, c in coefficients(f)
        !isone(numerator(c)) && push!(CF_set, numerator(c))
        !isone(denominator(c)) && push!(CF_set, denominator(c))
    end

    CF = sort!(collect(CF_set), rev=true)
    Lprim = ZZRingElem[]

    while length(Lprim) < N
        cur_prim = next_prime(rand(low:up))
        is_lucky = !(cur_prim in Lprim)
        idx = firstindex(CF)
        # Exploit decreasing order of CF
        while is_lucky && idx <= lastindex(CF) && CF[idx] > cur_prim
            is_lucky = !is_divisible_by(CF[idx], cur_prim)
            idx += 1
        end
        is_lucky && push!(Lprim, cur_prim)
    end
    return Lprim
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
        @assert ri < li1 "Intervals are not disjoint.", i, ri, li1
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

# Compute the image of a hyperbox B though a linear map L
# the intervals in B are supposed to be ordered as the variables in parent(L)
function _linear_interval_eval(L::QQMPolyRingElem,
                       B::Vector{Vector{QQFieldElem}})
    R = parent(L)
    x = gens(R)

    lo = hi = zero(QQ)

    for (xi, (a, b)) in zip(x, B)
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