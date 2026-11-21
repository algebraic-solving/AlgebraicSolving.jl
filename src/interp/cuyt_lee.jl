mutable struct Atomic{T}
    @atomic x::T
end

struct CuytLeeError <: Exception
    msg::String
end

Base.show(io::IO, e::CuytLeeError) = print(io, "CuytLeeError: ", e.msg)

@doc Markdown.doc"""
    _random_point(n::Int)

Generates a random point in $\mathbb{Z}^n$ with coordinates between 1 and 99.

**Note**: This is an internal function.
"""
function _random_point(n::Int)::Vector{Int}
    map(a -> 1 + abs(a) % 99, rand(Int, n))
end

@doc Markdown.doc"""
    _estimate_total_degree(R::QQMPolyRing, bb::Function; samples=5, show_progress=false)

Estimate the total degree of the rational function corresponding to the black box function `bb` by evaluating it at random points and performing univariate Thiele interpolation.

**Note**: This is an internal function.
"""
function _estimate_total_degree(
    R::QQMPolyRing,
    bb::Function;
    samples::Int=5,
    show_progress::Bool=false
)::Tuple{Int,Int}
    t = length(gens(R))
    @assert t > 0
    R_z, _ = polynomial_ring(QQ, :z)
    total_degree_counts = Dict{Tuple{Int,Int},Int}()
    total = 0
    while total < samples
        x = _random_point(t)
        try
            f = thiele(R_z, k -> bb(k .* x); show_progress=show_progress, offset=1)
            total_degree = (degree(denominator(f)), degree(numerator(f)))
            total_degree_counts[total_degree] = get(total_degree_counts, total_degree, 0) + 1
            total += 1
            if findmax(total_degree_counts)[1] > samples ÷ 2
                break
            end
        catch e
            if isa(e, ThieleError)
                continue
            else
                rethrow(e)
            end
        end
    end
    return findmax(total_degree_counts)[2]
end

@doc Markdown.doc"""
    _homogenize(f::QQMPolyRingElem, d::Int)

Homogenize the given polynomial `f` to total degree `d` by multiplying each term by the appropriate power of the first variable.

**Note**: This is an internal function.
"""
function _homogenize(
    f::QQMPolyRingElem,
    d::Int
)::QQMPolyRingElem
    R = parent(f)
    C = MPolyBuildCtx(R)
    for i in 1:f.length
        exp = exponent_vector(f, i)
        @assert exp[1] == 0
        total_deg = sum(exp)
        exp[1] += d - total_deg
        @assert exp[1] >= 0
        push_term!(C, coeff(f, i), exp)
    end
    return finish(C)
end

@doc Markdown.doc"""
    cuyt_lee_shifted(R::QQMPolyRing, bb::Function; <keyword arguments>)

Compute the multivariate rational function corresponding to the black box function `bb` using Cuyt and Lee's interpolation algorithm.
This function assumes that a random shift has already been applied to the input of the black box function, and does not perform any retries if the interpolation fails.

**Note**: This is an internal function. For a user-facing function that automatically applies random shifts, see `cuyt_lee`.
"""
function cuyt_lee_shifted(
    R::QQMPolyRing,
    bb::Function;
    retry::Int=10,
    nr_thrds::Int=1,
    show_progress::Bool=false,
    desc::String="Multivariate rational interpolation"
)::FracFieldElem{QQMPolyRingElem}
    # https://arxiv.org/pdf/1608.01902
    t = length(gens(R))
    if t == 0
        return R(bb(Vector{QQFieldElem}())) // one(R)
    end
    R_z, _ = polynomial_ring(QQ, :z)
    d_den, d_num = _estimate_total_degree(R, bb; show_progress=show_progress)
    D = max(d_den, d_num)
    m = t - 1
    idx = _simplex_indices(m, D)
    N = length(idx)
    nodes = [ZZ.(1:D+1) for _ in 1:m]
    point(a) = [ZZ(1); _simplex_point(nodes, a)]
    used = Vector{Union{Nothing,Vector{ZZRingElem}}}(nothing, N)
    coeffs_den = Vector{Vector{QQFieldElem}}(undef, N)
    coeffs_num = Vector{Vector{QQFieldElem}}(undef, N)
    done = falses(N)
    prog = ProgressBar(total=N; desc=desc, enabled=show_progress)
    update!(prog, 0)
    prog_lock = ReentrantLock()
    # Returns whether the evaluation at the k-th lattice point succeeded.
    function evaluate_point(k::Int, offset::Int)::Bool
        cur = point(idx[k])
        used[k] = cur
        f = nothing
        try
            f = thiele(R_z, z -> bb(z .* cur); retry=retry, show_progress=show_progress, offset=offset)
        catch e
            if isa(e, BoundsError) || isa(e, ThieleError)
                return false
            else
                rethrow(e)
            end
        end
        if degree(denominator(f)) != d_den || degree(numerator(f)) != d_num
            return false
        end
        c = constant_coefficient(denominator(f))
        if c == 0
            return false
        end
        coeffs_den[k] = collect(coefficients(denominator(f))) ./ c
        coeffs_num[k] = collect(coefficients(numerator(f))) ./ c
        return true
    end
    # A failing point is unlucky for one of its coordinates,
    # so we replace one of its nodes and recompute the lattice points sharing that node.
    max_repairs = retry * (D + 1)
    repairs = 0
    while true
        pending = findall(!, done)
        isempty(pending) && break
        ok = falses(length(pending))
        nt = max(1, min(nr_thrds, length(pending)))
        Threads.@threads for c in 1:nt
            for l in c:nt:length(pending)
                ok[l] = evaluate_point(pending[l], c)
                if ok[l]
                    lock(() -> next!(prog), prog_lock)
                end
            end
        end
        done[pending[ok]] .= true
        for l in findall(!, ok)
            k = pending[l]
            a = idx[k]
            # skip points whose coordinates were already changed by an earlier repair
            used[k] == point(a) || continue
            repairs += 1
            # with a single variable there is no node to change: retry with another shift
            if m == 0 || repairs > max_repairs
                finish!(prog)
                throw(CuytLeeError("Failed to collect enough data points for interpolation. This could happen if the black box function is singular at zero, or if the expected total degree is incorrect."))
            end
            j, aj = _replace_node!(nodes, a)
            for k′ in 1:N
                if idx[k′][j] == aj && done[k′]
                    done[k′] = false
                end
            end
        end
        update!(prog, count(done))
    end
    # We interpolate the homogeneous parts of the denominator and numerator separately,
    # each on the sub-lattice of the matching total degree.
    den = zero(R)
    num = zero(R)
    for i in 0:d_den
        y = [coeffs_den[k][i+1] for k in 1:binomial(m + i, m)]
        den += _homogenize(_newton_simplex(R, collect(2:t), nodes, y, i), i)
    end
    for i in 0:d_num
        y = [coeffs_num[k][i+1] for k in 1:binomial(m + i, m)]
        num += _homogenize(_newton_simplex(R, collect(2:t), nodes, y, i), i)
    end
    finish!(prog)
    return num // den
end

@doc Markdown.doc"""
    cuyt_lee_with_shift(R::QQMPolyRing, bb::Function, shift::Vector{Int}; <keyword arguments>)

Compute the multivariate rational function corresponding to the black box function `bb` using Cuyt and Lee's interpolation algorithm,
with a given shift applied to the input of the black box function.

**Note**: This is an internal function. For a user-facing function that automatically applies random shifts, see `cuyt_lee`.
"""
function cuyt_lee_with_shift(
    R::QQMPolyRing,
    bb::Function,
    shift::Vector{Int};
    retry::Int=10,
    nr_thrds::Int=1,
    show_progress::Bool=false,
    desc::String="Multivariate rational interpolation"
)::FracFieldElem{QQMPolyRingElem}
    t = length(gens(R))
    if t == 0
        return R(bb(Vector{QQFieldElem}())) // one(R)
    end
    f_shifted = cuyt_lee_shifted(R, z -> bb(z .+ shift); retry=retry, nr_thrds=nr_thrds, show_progress=show_progress, desc=desc)
    x = gens(R) .- shift
    num = evaluate(numerator(f_shifted), x)
    den = evaluate(denominator(f_shifted), x)
    return num // den
end

@doc Markdown.doc"""
    cuyt_lee(R::QQMPolyRing, bb::Function; <keyword arguments>)

Compute the multivariate rational function corresponding to the black box function `bb` using Cuyt and Lee's interpolation algorithm.

# Arguments
- `R::QQMPolyRing`: the multivariate polynomial ring over the rationals.
- `bb::Function`: a black box function that takes a vector of rational numbers as input and returns a rational number as output.
- `initial_shift::Vector{Int}=_random_point(length(gens(R)))`: the initial shift to use for the interpolation.
- `retry::Int=10`: the maximum number of consecutive failures allowed when evaluating the black box function or interpolating the points.
- `nr_thrds::Int=1`: the number of threads to use when evaluating the black box function.
- `show_progress::Bool=false`: whether to show a progress bar while collecting points.
- `desc::String="Multivariate rational interpolation"`: the description to show in the progress bar.

# Example
```jldoctest
julia> using AlgebraicSolving

julia> using AlgebraicSolving.Interpolation

julia> R, _ = polynomial_ring(QQ, ["a", "b"]);

julia> f = t -> (t[1] + t[2])^2 / (t[1] - t[2])^2;

julia> cuyt_lee(R, f)
(a^2 + 2*a*b + b^2)//(a^2 - 2*a*b + b^2)
```
"""
function cuyt_lee(
    R::QQMPolyRing,
    bb::Function;
    initial_shift=_random_point(length(gens(R))),
    retry::Int=10,
    nr_thrds::Int=1,
    show_progress::Bool=false,
    desc::String="Multivariate rational interpolation"
)::FracFieldElem{QQMPolyRingElem}
    t = length(gens(R))
    if t == 0
        return R(bb(Vector{QQFieldElem}())) // one(R)
    end
    shift = initial_shift
    for i in 1:retry
        try
            return cuyt_lee_with_shift(R, bb, shift; retry=retry, nr_thrds=nr_thrds, show_progress=show_progress, desc=desc)
        catch e
            if isa(e, CuytLeeError)
                if show_progress
                    @warn "Interpolation failed, retrying with a different shift... Retries left: $(retry - i)"
                end
                shift = _random_point(t)
            else
                rethrow(e)
            end
        end
    end
    throw(CuytLeeError("Interpolation failed after maximum number of retries."))
end
