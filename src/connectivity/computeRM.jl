#using Pkg
#Pkg.activate("AlgebraicSolving.jl")
#using Revise
#using AlgebraicSolving
#using Nemo

export roadmap, computepolar, MidRationalPoints, fbr, all_eqs, all_base_pts, nb_nodes, compute_minors, compute_minors_bis, computepolarL
include("Cannytools.jl")

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
- `q::Vector{QQFieldElem}=QQFieldElem[]`: single base point with rational coefficients
- `C::Vector{Vector{QQFieldElem}}=Vector{QQFieldElem}[]`: query points with rational coefficients
- `info_level::Int=0`: verbosity level
- `checks::Bool=false`: whether perform checks (dimension, regularity, etc.)
- 'generic_change=false": whether it performs a prior random linear change of variables (TODO)
)

# Examples
```jldoctest
julia> using AlgebraicSolving

julia> R,(x1,x2,x3,x4) = polynomial_ring(QQ, ["x1","x2","x3","x4"])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x1, x2, x3, x4])

julia> I = Ideal([x1^2+x2^2+x3^2+x4^2+9-1)^2-4*9*(x1^2+x2^2+x3^2)])
QQMPolyRingElem[x1^4 + 2*x1^2*x2^2 + 2*x1^2*x3^2 + 2*x1^2*x4^2 - 20*x1^2 + x2^4 + 2*x2^2*x3^2 + 2*x2^2*x4^2 - 20*x2^2 + x3^4 + 2*x3^2*x4^2 - 20*x3^2 + x4^4 + 16*x4^2 + 64]

julia> RM = roadmap(I)
Vector{QQFieldElem}[[], [-3], [-3, -2], [-1], [-1, -2], [-1, -1], [-1, 2], [3], [3, 2]]

julia> nb_nodes(RM)
9

julia> all_eqs(RM)
9-element Vector{Ideal{QQMPolyRingElem}}:
 QQMPolyRingElem[x1^4 + 2*x1^2*x2^2 + 2*x1^2*x3^2 + 2*x1^2*x4^2 - 20*x1^2 + x2^4 + 2*x2^2*x3^2 + 2*x2^2*x4^2 - 20*x2^2 + x3^4 + 2*x3^2*x4^2 - 20*x3^2 + x4^4 + 16*x4^2 + 64, 4*x1^2*x3 + 4*x2^2*x3 + 4*x3^3 + 4*x3*x4^2 - 40*x3, 4*x1^2*x4 + 4*x2^2*x4 + 4*x3^2*x4 + 4*x4^3 + 32*x4]
 QQMPolyRingElem[x1^4 + 2*x1^2*x2^2 + 2*x1^2*x3^2 + 2*x1^2*x4^2 - 20*x1^2 + x2^4 + 2*x2^2*x3^2 + 2*x2^2*x4^2 - 20*x2^2 + x3^4 + 2*x3^2*x4^2 - 20*x3^2 + x4^4 + 16*x4^2 + 64, 4*x1^2*x4 + 4*x2^2*x4 + 4*x3^2*x4 + 4*x4^3 + 32*x4, x1 + 3]
 QQMPolyRingElem[x1^4 + 2*x1^2*x2^2 + 2*x1^2*x3^2 + 2*x1^2*x4^2 - 20*x1^2 + x2^4 + 2*x2^2*x3^2 + 2*x2^2*x4^2 - 20*x2^2 + x3^4 + 2*x3^2*x4^2 - 20*x3^2 + x4^4 + 16*x4^2 + 64, x1 + 3, x2 + 2]
 QQMPolyRingElem[x1^4 + 2*x1^2*x2^2 + 2*x1^2*x3^2 + 2*x1^2*x4^2 - 20*x1^2 + x2^4 + 2*x2^2*x3^2 + 2*x2^2*x4^2 - 20*x2^2 + x3^4 + 2*x3^2*x4^2 - 20*x3^2 + x4^4 + 16*x4^2 + 64, 4*x1^2*x4 + 4*x2^2*x4 + 4*x3^2*x4 + 4*x4^3 + 32*x4, x1 + 1]
 QQMPolyRingElem[x1^4 + 2*x1^2*x2^2 + 2*x1^2*x3^2 + 2*x1^2*x4^2 - 20*x1^2 + x2^4 + 2*x2^2*x3^2 + 2*x2^2*x4^2 - 20*x2^2 + x3^4 + 2*x3^2*x4^2 - 20*x3^2 + x4^4 + 16*x4^2 + 64, x1 + 1, x2 + 2]
 QQMPolyRingElem[x1^4 + 2*x1^2*x2^2 + 2*x1^2*x3^2 + 2*x1^2*x4^2 - 20*x1^2 + x2^4 + 2*x2^2*x3^2 + 2*x2^2*x4^2 - 20*x2^2 + x3^4 + 2*x3^2*x4^2 - 20*x3^2 + x4^4 + 16*x4^2 + 64, x1 + 1, x2 + 1]
 QQMPolyRingElem[x1^4 + 2*x1^2*x2^2 + 2*x1^2*x3^2 + 2*x1^2*x4^2 - 20*x1^2 + x2^4 + 2*x2^2*x3^2 + 2*x2^2*x4^2 - 20*x2^2 + x3^4 + 2*x3^2*x4^2 - 20*x3^2 + x4^4 + 16*x4^2 + 64, x1 + 1, x2 - 2]
 QQMPolyRingElem[x1^4 + 2*x1^2*x2^2 + 2*x1^2*x3^2 + 2*x1^2*x4^2 - 20*x1^2 + x2^4 + 2*x2^2*x3^2 + 2*x2^2*x4^2 - 20*x2^2 + x3^4 + 2*x3^2*x4^2 - 20*x3^2 + x4^4 + 16*x4^2 + 64, 4*x1^2*x4 + 4*x2^2*x4 + 4*x3^2*x4 + 4*x4^3 + 32*x4, x1 - 3]
 QQMPolyRingElem[x1^4 + 2*x1^2*x2^2 + 2*x1^2*x3^2 + 2*x1^2*x4^2 - 20*x1^2 + x2^4 + 2*x2^2*x3^2 + 2*x2^2*x4^2 - 20*x2^2 + x3^4 + 2*x3^2*x4^2 - 20*x3^2 + x4^4 + 16*x4^2 + 64, x1 - 3, x2 - 2]
```
"""
function roadmap(
    I::Ideal{T} where T <: QQMPolyRingElem;                 # input ideal
    q::Vector{QQFieldElem}=QQFieldElem[],                   # single base point with rational coefficients
    C::Vector{Vector{QQFieldElem}}=Vector{QQFieldElem}[],   # query points with rational coefficients
    info_level::Int=0,                                               # verbosity level
    checks::Bool=false                                      # perform checks (dimension, regularity, etc.)
)
    # Some base cases
    if nvars(parent(I))<=2
        return [I.gens]
    end
    # Some preprocessing
    isnothing(I.dim) && dimension(I)
    # Base points
    e = length(q)

    ## Fq ##
    # Genericity assumption (can be checked)
    if checks
        @assert(dimension(fbr(I,q)) == I.dim - e, "Non-generic coordinates")
    end

    # Terminal case (dim <=1)
    if I.dim - e <= 1
        return RMnode([], q, RMnode[])
    end

    ## sing(Fq) ##
    if checks
        info_level>0 && println("Check real quasi-smoothness")
        singFq = computepolar(1:e, I) |> Ideal
        @assert(isempty(real_solutions(fbr(singFq, q), info_level=max(info_level-1,0), nr_thrds=Threads.nthreads())),
                "Non-empty real sing locus!")
    end

    ## K(pi_1,Fq) ##
    info_level>0 && println("I-critical points")
    K1Fq = computepolar(1:e+1, I) |> Ideal
    K1Fq = real_solutions(fbr(K1Fq,q), info_level=max(info_level-1,0), nr_thrds=Threads.nthreads(), interval=true)

    ## K(pi_2, Fq) ##
    info_level>0 && println("Polar variety")
    K2Fqmins = computepolar(1:e+2, I, only_mins=true)
    K2Fq = vcat(I.gens, K2Fqmins) |> Ideal
    if checks
        @assert(fbr(K2Fq, q) |> dimension == 1, "Non-generic polar variety")
    else
        K2Fq.dim = e + 1
    end
    RM = RMnode(K2Fqmins, q, RMnode[])

    ## Points with vertical tg in K(pi_2, Fq) ##
    info_level>0 && println("W-critical points with vertical tangent")
    K1WmFq = computepolar(1:e+2, K2Fq, dimproj=e) |> Ideal
    K1WmFq = real_solutions(fbr(K1WmFq,q), info_level=max(info_level-1,0), nr_thrds=Threads.nthreads(), interval=true)

    ## New base and query points ##
    Cq = isempty(q) ? C : [ c for c in C if c[e] == q[e]]
    K1W = vcat(K1Fq, K1WmFq)
    # Heuristic to be proven (Reeb's th)
    #K1W = K1W[2:end-1]
    ##########
    K1WRat = MidRationalPoints(getindex.(K1W,e+1), unique(getindex.(Cq, e+1)))
    newQ = vcat.(Ref(q), K1WRat)

    # Recursively compute roadmap of possible fibers
    if !isempty(newQ)
        for newq in newQ
            RMFq = roadmap(I, q=newq, C=Cq)
            push!(RM.children, RMFq)
        end
    end

    if e == 0
        return Roadmap(I, RM)
    else
        return RM
    end
end

function roadmap(
    I::Ideal{P},                # input ideal
    C::Ideal{P};                # ideal defining query points
    info_level::Int=0,                   # verbosity level
    checks::Bool=false          # perform checks (dimension, regularity, etc.)
) where (P <: QQMPolyRingElem)
    @assert(parent(I)==parent(C), "Equations for variety and query points must live the same ring")
    CQ = real_solutions(C, info_level=max(info_level-1,0), nr_thrds=Threads.nthreads())
    return roadmap(I, C=CQ, info_level=info_level, checks=checks)
end


function fbr(F::Vector{P} where P <: QQMPolyRingElem, Q::Vector{QQFieldElem})
    @assert(!isempty(F), "Empty polynomial vector")
    vars = gens(parent(first(F)))
    return Ideal(vcat(F, [vars[i] - Q[i] for i in 1:min(length(vars),length(Q))]))
end

function fbr(I::Ideal{P} where P <: QQMPolyRingElem, Q::Vector{QQFieldElem})
    return fbr(I.gens, Q)
end

#=
## Test ##
R,(x1,x2,x3,x4) = polynomial_ring(QQ, ["x1","x2","x3","x4"])
#h = x1^2+x2^2+x3^2+x4^2-1
h = (x1^2+x2^2+x3^2+x4^2+9-1)^2-4*9*(x1^2+x2^2+x3^2)
h = evaluate(h,[x1,x2,x3],[x1+rand(-10:10)*x2+rand(-10:10)*x3+rand(-10:10)*x4,x2+rand(-10:10)*x3+rand(-10:10)*x4,x3+rand(-10:10)*x4])
I = Ideal([h])

carte = (computeRM(I, 3, [Vector{QQFieldElem}(undef,0)]))
#println(carte)
=#


