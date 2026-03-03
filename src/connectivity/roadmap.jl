#using Pkg
#Pkg.activate("AlgebraicSolving.jl")
#using Revise
#using AlgebraicSolving
#using Nemo

export roadmap, computepolar, _mid_rational_points, _mid_rational_points_inter, _fbr, all_eqs, all_base_pts, nb_nodes, computepolar
include("polar.jl")

@doc Markdown.doc"""
    roadmap(I::Ideal{T} where T <: MPolyRingElem, <keyword arguments>)

Given a **radical** ideal `I` with solution set X, that is smooth and
whose real trace XR is bounded, return a roadmap of XR

The output is given as a Roadmap structure, encoding the recursive structure
of roadmaps. It is encoded as a chained list, whose root is containing the equations defining X
and each node representing a curve component, that is defined by additional polar equation and base point<
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
) where (P <: QQMPolyRingElem)
    return _roadmap_rec(I, QQFieldElem[], C, info_level, checks)
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
) where (P <: QQMPolyRingElem)
    @assert(parent(I)==parent(C), "Equations for variety and query points must live the same ring")
    CI = real_solutions(C, info_level=max(info_level-1,0), nr_thrds=Threads.nthreads())
    return _roadmap_rec(I, QFieldElem[], CI, info_level, checks)
end

function _roadmap_rec(
    I::Ideal{T} where T <: QQMPolyRingElem,     # input ideal
    q::Vector{QQFieldElem},                     # single base point with rational coefficients
    C::Vector{Vector{Vector{QQFieldElem}}},     # query points with rational coefficients
    info_level::Int,                            # verbosity level
    checks::Bool                                # perform checks (dimension, regularity, etc.)
)
    # Some base cases
    if nvars(parent(I))<=2
        return [I.gens]
    end

    # Some preprocessng
    if isnothing(I.dim)
        lucky_prime = _generate_lucky_primes(I.gens, one(ZZ)<<30, one(ZZ)<<31-1, 1) |> first
        local INEW = Ideal(change_base_ring.(Ref(GF(lucky_prime)), I.gens))
        dimension(INEW)
        I.dim = INEW.dim
    end
    e = length(q)

    ## Fq ##
    # Genericity assumption (can be checked)
    if checks
        local Fnew = _fbr(I,q).gens
        new_lucky_prime = _generate_lucky_primes(Fnew, one(ZZ)<<30, one(ZZ)<<31-1, 1) |> first
        local INEW = Ideal(change_base_ring.(Ref(GF(new_lucky_prime)), Fnew))
        @assert(dimension(INEW) == I.dim - e, "Non-generic coordinates")
    end

    # Terminal case (dim <=1)
    if I.dim - e <= 1
        return RMnode(q, [], RMnode[])
    end

    ## sing(Fq) ##
    if checks
        info_level>0 && println("Check smoothness")
        local FNEW = _fbr(computepolar(1:e, I)|> Ideal, q).gens
        new_lucky_prime = _generate_lucky_primes(FNEW, one(ZZ)<<30, one(ZZ)<<31-1, 1) |> first
        local INEW = Ideal(change_base_ring.(Ref(GF(new_lucky_prime)), FNEW))
        @assert(dimension(INEW) == -1, "Non-empty sing locus!")
    end

    ## K(pi_1,Fq) ##
    info_level>0 && println("Compute x1-critical points: K1")
    K1Fq = computepolar(1:e+1, I) |> Ideal
    K1Fq = real_solutions(_fbr(K1Fq,q), info_level=max(info_level-1,0), nr_thrds=Threads.nthreads(), interval=true)

    ## K(pi_2, Fq) ##
    info_level>0 && println("Compute (x1,x2)-polar variety: W")
    K2Fqmins = computepolar(1:e+2, I, only_mins=true)
    K2Fq = vcat(I.gens, K2Fqmins) |> Ideal
    if checks
        local FNEW = _fbr(K2Fq, q).gens
        new_lucky_prime = _generate_lucky_primes(FNEW, one(ZZ)<<30, one(ZZ)<<31-1, 1) |> first
        local INEW = Ideal(change_base_ring.(Ref(GF(new_lucky_prime)), FNEW))
        @assert(dimension(INEW) == 1, "Non-generic polar variety")
    else
        K2Fq.dim = e + 1
    end
    RM = RMnode(q, K2Fqmins, RMnode[])

    ## Points with vertical tg in K(pi_2, Fq) ##
    info_level>0 && println("Compute W-critical points with vertical tangent: K1W")
    K1WmFq = computepolar(1:e+2, K2Fq, dimproj=e) |> Ideal
    K1WmFq = real_solutions(_fbr(K1WmFq,q), info_level=max(info_level-1,0), nr_thrds=Threads.nthreads(), interval=true)

    ## New base and query points ##
    Cq = isempty(q) ? C : [ c for c in C if q[e] in c[e]]
    K1W = vcat(K1Fq, K1WmFq)
    # Heuristic to be proven (Reeb's th)
    #K1W = K1W[2:end-1]
    ##########
    K1WRat = _mid_rational_points_inter(getindex.(K1W,e+1), unique(getindex.(Cq, e+1)))
    newQ = vcat.(Ref(q), K1WRat)

    # Recursively compute roadmap of possible fibers
    if !isempty(newQ)
        for newq in newQ
            RMFq = _roadmap_rec(I, newq, Cq, info_level, checks)
            push!(RM.children, RMFq)
        end
    end

    if e == 0
        return Roadmap(I, RM)
    else
        return RM
    end
end

function _mid_rational_points_inter(S::Vector{Vector{T}}, Q::Vector{Vector{T}} = Vector{T}[]) where {T <: QQFieldElem}
    # * S is a list of [ [l_1,r_1], ..., [l_n, r_n] ]
    # such that the [l_i, r_i] are rational and disjoint open intervals.
    # * Same assumptions on Q
    # * Intervals in S and Q do not intersect as well
    #
    # It orders the [l_i,r_i], and compute a list ratioP
    # intersecting all intervals of Q and such that
    # strictly between each of the [l_i,r_i] there is:
    # - either at least one element inside an interval of Q
    # - or the simplest rational number
    # TODO:
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


