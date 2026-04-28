@doc Markdown.doc"""
    _to_multivariate(R::QQMPolyRing, p::QQPolyRingElem)::QQMPolyRingElem

Convert a univariate polynomial `p` to a multivariate polynomial in the ring `R` with the same variable.

**Note**: This is an internal function.
"""
function _to_multivariate(R::QQMPolyRing, p::QQPolyRingElem)::QQMPolyRingElem
    @assert length(gens(R)) == 1
    M = MPolyBuildCtx(R)
    for i in 0:length(p)
        push_term!(M, coeff(p, i), [i])
    end
    finish(M)
end

@doc Markdown.doc"""
    _sample_points(f::Vector{QQPolyRingElem}; worker_pool::AbstractWorkerPool=default_worker_pool())::Vector{QQFieldElem}

Given a vector of univariate polynomials, returns a vector of rational numbers such that each connected component of the non-vanishing set of the polynomials contains at least one of the returned points.

**Note**: This is an internal function.
"""
function _sample_points(f::Vector{QQPolyRingElem}; worker_pool::AbstractWorkerPool=default_worker_pool())::Vector{QQFieldElem}
    R, _ = polynomial_ring(QQ, [:x])
    fₓ = map(p -> _to_multivariate(R, p), f)
    factors = unique([p for fₓ′ in fₓ for (p, _) in factor(fₓ′)])
    # We map each factor to its roots, and each root to its factor
    roots_by_factor = Dict{QQMPolyRingElem,Vector{Vector{Vector{QQFieldElem}}}}()
    factors_by_root = Dict{Vector{Vector{QQFieldElem}},QQMPolyRingElem}()
    for factor in factors
        roots = real_solutions(Ideal([factor]); interval=true, worker_pool=worker_pool)
        roots_by_factor[factor] = roots
        for r in roots
            factors_by_root[r] = factor
        end
    end
    # We order intervals by their left endpoint
    roots = sort(collect(keys(factors_by_root)), by=x -> x[1][1])
    # If intervals are not ordered by their right endpoint,
    # we merge the offending factors
    while (!issorted(roots, by=x -> x[1][2]))
        i = findfirst(i -> roots[i][1][2] > roots[i+1][1][2], 1:length(roots)-1)
        bad_factors = [factors_by_root[roots[i]], factors_by_root[roots[i+1]]]
        # Remove old factors and roots
        for bad_factor in bad_factors
            for r in roots_by_factor[bad_factor]
                delete!(factors_by_root, r)
            end
            delete!(roots_by_factor, bad_factor)
        end
        # Replace with merged factor and its roots
        merged_factor = prod(bad_factors)
        merged_roots = real_solutions(Ideal([merged_factor]); interval=true, worker_pool=worker_pool)
        roots_by_factor[merged_factor] = merged_roots
        for r in merged_roots
            factors_by_root[r] = merged_factor
        end
        roots = sort(collect(keys(factors_by_root)), by=r -> r[1][1])
    end
    # Now the intervals are properly ordered, we can sample points
    points = Vector{QQFieldElem}()
    if length(roots) == 0
        push!(points, QQ(0))
        return points
    end
    push!(points, floor(roots[1][1][1]) - 1)
    for i in 1:length(roots)-1
        push!(points, (roots[i][1][2] + roots[i+1][1][1]) // 2)
    end
    push!(points, ceil(roots[end][1][2]) + 1)
    return points
end

function _sample_points_0(f::Vector{QQMPolyRingElem})::Vector{Vector{QQFieldElem}}
    @assert all(map(is_constant, f))
    return [Vector{QQFieldElem}()]
end

function _sample_points_1(f::Vector{QQMPolyRingElem}; worker_pool::AbstractWorkerPool=default_worker_pool())::Vector{Vector{QQFieldElem}}
    @assert all(map(is_univariate, f))
    R, _ = polynomial_ring(QQ, :x)
    return [[p] for p in _sample_points(map(p -> to_univariate(R, p), f); worker_pool=worker_pool)]
end

function _sample_points_desc(n::Int)::String
    if n == 0
        return "Constant"
    elseif n == 1
        return "Univariate"
    elseif n == 2
        return "Bivariate"
    elseif n == 3
        return "Trivariate"
    else
        return "Multivariate"
    end
end

function _sample_points_2(
    f::Vector{QQMPolyRingElem},
    xs::Vector{QQMPolyRingElem};
    nr_thrds::Int=1,
    worker_pool::AbstractWorkerPool=default_worker_pool(),
    show_progress::Bool=false,
    desc::String="$(_sample_points_desc(length(xs))) sample points"
)::Vector{Vector{QQFieldElem}}
    @assert length(xs) >= 2
    x₁ = xs[1:end-1]
    x₂ = xs[end]
    factors = unique([p for fₓ in f for (p, _) in factor(fₓ)])
    v = Vector{QQMPolyRingElem}()
    for i in eachindex(factors)
        if !isempty(intersect(x₁, vars(factors[i])))
            push!(v, leading_coefficient(factors[i], length(xs)))
        end
        if x₂ in vars(factors[i])
            push!(v, Interpolation.discriminant(factors[i], length(xs); nr_thrds=nr_thrds))
        end
    end
    for i in 1:length(factors)-1
        for j in i+1:length(factors)
            if x₂ in vars(factors[i]) || x₂ in vars(factors[j])
                push!(v, Interpolation.resultant(factors[i], factors[j], length(xs); nr_thrds=nr_thrds))
            end
        end
    end
    p₁ = if length(xs) == 2
        _sample_points_1(v; worker_pool=worker_pool)
    else
        _sample_points_2(v, x₁; show_progress=show_progress, worker_pool=worker_pool)
    end
    prog = Progress.ProgressBar(total=length(p₁); desc=desc, enabled=show_progress)
    Progress.update!(prog, 0)
    function _points_chunk(p::Vector{Vector{QQFieldElem}})::Vector{Vector{QQFieldElem}}
        res_chunk = Vector{Vector{QQFieldElem}}()
        for i in eachindex(p)
            p₂ = _sample_points_1(map(p′ -> evaluate(p′, x₁, p[i]), f); worker_pool=worker_pool)
            for j in eachindex(p₂)
                push!(res_chunk, [p[i]; p₂[j][1]])
            end
            Progress.next!(prog)
        end
        res_chunk
    end
    chunk_size = ceil(Int, length(p₁) / nr_thrds)
    chunks = [p₁[i:min(i + chunk_size - 1, end)] for i in 1:chunk_size:length(p₁)]
    tasks = [Threads.@spawn _points_chunk(chunk) for chunk in chunks]
    res = sort(vcat(fetch.(tasks)...))
    Progress.finish!(prog)
    res
end

@doc Markdown.doc"""
    _sample_points(f::Vector{QQMPolyRingElem}; nr_thrds::Int=1, worker_pool::AbstractWorkerPool=default_worker_pool(), show_progress::Bool=false)::Vector{Vector{QQFieldElem}}

Given a vector of multivariate polynomials, returns a vector of points such that each connected component of the non-vanishing set of the polynomials contains at least one of the returned points.

**Note**: This is an internal function.
"""
function _sample_points(
    f::Vector{QQMPolyRingElem};
    nr_thrds::Int=1,
    worker_pool::AbstractWorkerPool=default_worker_pool(),
    show_progress::Bool=false
)::Vector{Vector{QQFieldElem}}
    if any(is_zero, f)
        error("Cannot sample points for polynomials containing zero polynomial")
    end
    p = Vector{Vector{QQFieldElem}}([])
    if length(gens(parent(f[1]))) == 0
        p = _sample_points_0(f)
    elseif length(gens(parent(f[1]))) == 1
        p = _sample_points_1(f; worker_pool=worker_pool)
    else
        p = _sample_points_2(f, gens(parent(f[1])); nr_thrds=nr_thrds, worker_pool=worker_pool, show_progress=show_progress)
    end
    return p
end

@doc Markdown.doc"""
    points_per_components(eqs::Vector{QQMPolyRingElem}, pos::Vector{QQMPolyRingElem}, ineqs::Vector{QQMPolyRingElem}; nr_thrds::Int=1, worker_pool::AbstractWorkerPool=default_worker_pool(), info_level::Int=0)::Vector{Vector{Vector{QQFieldElem}}}

Given a semi-algebraic set defined by equations `eqs`, positivity constraints `pos`, and non-vanishing constraints `ineqs`, returns a vector of points meeting all connected components of the set.
Each point is represented by an isolating box, i.e., a vector of intervals represented by pairs of rational numbers.

**NOTE**: Sampling points per components is currently only implemented for systems without equations, i.e., `eqs` must be empty.

# Arguments
- `eqs::Vector{QQMPolyRingElem}`: equations defining the semi-algebraic set.
- `pos::Vector{QQMPolyRingElem}`: positivity constraints defining the semi-algebraic set.
- `ineqs::Vector{QQMPolyRingElem}`: non-vanishing constraints defining the semi-algebraic set.
- `nr_thrds::Int=1`: the number of threads to use for parallel computations.
- `worker_pool::AbstractWorkerPool=default_worker_pool()`: the worker pool to use for parallel computations.
- `info_level::Int=0`: info level printout: off (`0`, default), summary (`1`), detailed (`2`).

# Examples

```jldoctest
julia> using AlgebraicSolving

julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> length(points_per_components(QQMPolyRingElem[], [x, y], [x + y - 1])) >= 2
true
```
"""
function points_per_components(
    eqs::Vector{QQMPolyRingElem},
    pos::Vector{QQMPolyRingElem},
    ineqs::Vector{QQMPolyRingElem};
    nr_thrds::Int=1,
    worker_pool::AbstractWorkerPool=default_worker_pool(),
    info_level::Int=0
)::Vector{Vector{Vector{QQFieldElem}}}
    if !isempty(eqs)
        throw(ArgumentError("Sampling points per components is not implemented for systems involving equations"))
    end
    points = _sample_points([pos; ineqs]; nr_thrds=nr_thrds, worker_pool=worker_pool, show_progress=info_level > 0)
    res = Vector{Vector{Vector{QQFieldElem}}}()
    for point in points
        if all(evaluate(constraint, point) > 0 for constraint in pos)
            push!(res, map(x -> [x, x], point))
        end
    end
    return res
end
