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
    d = [0; fill(max(d_den, d_num), t - 1)...]
    x = Vector{Vector{ZZRingElem}}()
    coeffs_den = []
    coeffs_num = []
    data_lock = ReentrantLock()
    prog = ProgressBar(total=prod(d .+ 1); desc=desc, enabled=show_progress)
    update!(prog, 0)
    function populate(cur::Vector{ZZRingElem}, dim::Int; num_threads::Int=1, offset=1)::Bool
        if dim > t
            f = nothing
            try
                f = thiele(R_z, k -> bb(k .* cur); retry=retry, show_progress=show_progress, offset=offset)
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
            lock(data_lock) do
                push!(x, copy(cur))
                push!(coeffs_den, collect(coefficients(denominator(f))) ./ c)
                push!(coeffs_num, collect(coefficients(numerator(f))) ./ c)
                update!(prog, length(x))
            end
            return true
        end
        total = Atomic(0)
        failures = Atomic(0)
        i = 1
        while total.x < d[dim] + 1
            num_threads_chunk = min(num_threads, d[dim] + 1 - total.x)
            Threads.@threads for j in 0:num_threads_chunk-1
                if populate([cur; ZZ(i + j)], dim + 1; num_threads=1, offset=max(offset, j + 1))
                    @atomic total.x += 1
                    @atomic failures.x = 0
                else
                    @atomic failures.x += 1
                end
            end
            i += num_threads_chunk
            if failures.x >= retry
                return false
            end
        end
        return true
    end
    res = populate([ZZ(1)], 2; num_threads=nr_thrds)
    if !res
        finish!(prog)
        throw(CuytLeeError("Failed to collect enough data points for interpolation. This could happen if the black box function is singular at zero, or if the expected total degree is incorrect."))
    end
    perm = sortperm(x)
    x = x[perm]
    coeffs_den = coeffs_den[perm]
    coeffs_num = coeffs_num[perm]
    # We interpolate the denominator and numerator separately
    den = zero(R)
    num = zero(R)
    for i in 0:d_den
        y = [coeffs_den[j][i+1] for j in 1:length(x)]
        den += _homogenize(newton(R, x, y, d), i)
    end
    for i in 0:d_num
        y = [coeffs_num[j][i+1] for j in 1:length(x)]
        num += _homogenize(newton(R, x, y, d), i)
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
