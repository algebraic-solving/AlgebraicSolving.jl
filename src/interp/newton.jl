struct Fmpq
    num::Int
    den::Int
end

@doc Markdown.doc"""
    _flint_lagrange(R::QQPolyRing, x::Vector{ZZRingElem}, y::Vector{QQFieldElem})

Compute the univariate polynomial corresponding to the given points `x` and values `y` using Flint's implementation of Lagrange interpolation.

**Note**: This is an internal function.
"""
function _flint_lagrange(
    R::QQPolyRing,
    x::Vector{ZZRingElem},
    y::Vector{QQFieldElem}
)::QQPolyRingElem
    @assert length(x) == length(y)
    n = length(x)
    z = zero(R)
    ax = Vector{Int}(map(i -> Int(x[i].d), 1:n))
    ay = Vector{Fmpq}(map(i -> Fmpq(Int(y[i].num), Int(y[i].den)), 1:n))
    @ccall Nemo.libflint.fmpq_poly_interpolate_fmpz_fmpq_vec(z::Ref{QQPolyRingElem}, ax::Ptr{Int}, ay::Ptr{Fmpq}, n::Int)::Nothing
    @assert evaluate(z, x[n]) == y[n]
    return z
end

@doc Markdown.doc"""
    newton(R::QQMPolyRing, x::Vector{Vector{ZZRingElem}}, y::Vector{QQFieldElem}, d::Vector{Int})

Compute the multivariate polynomial corresponding to the given points `x` and values `y` of degrees at most `d` using Newton's interpolation formula.

**Note**: This is an internal function.
"""
function newton(
    R::QQMPolyRing,
    x::Vector{Vector{ZZRingElem}},
    y::Vector{QQFieldElem},
    d::Vector{Int}
)::QQMPolyRingElem
    n = length(x)
    @assert length(y) == prod(d .+ 1) == n
    step = prod(d[2:end] .+ 1)
    x_index = length(gens(R)) - length(d) + 1
    z = gens(R)[x_index]
    x_trunc = [x[i][x_index] for i in 1:step:n]
    if length(d) == 1
        R′, _ = polynomial_ring(QQ, :z)
        y_trunc = [y[i] for i in 1:(d[1]+1)]
        f = _flint_lagrange(R′, x_trunc, y_trunc)
        C = MPolyBuildCtx(R)
        for i in 0:length(f)-1
            c = coeff(f, i)
            if c != 0
                exp = fill(0, length(gens(R)))
                exp[x_index] = i
                push_term!(C, c, exp)
            end
        end
        res = finish(C)
        return res
    end
    y_trunc = [newton(R, x[i:i+step-1], y[i:i+step-1], d[2:end]) for i in 1:step:n]
    if d[1] == 0
        return y_trunc[1]
    end
    dd = [copy(y_trunc)]
    for j in 1:d[1]
        dd_j = Vector{QQMPolyRingElem}()
        for i in 1:(d[1]-j+1)
            push!(dd_j, (dd[j][i+1] - dd[j][i]) / (R(x_trunc[i+j]) - R(x_trunc[i])))
        end
        push!(dd, dd_j)
    end
    f = zero(R)
    term = one(R)
    for j in 0:d[1]
        f += dd[j+1][1] * term
        term *= (z - x_trunc[j+1])
    end
    return f
end

@doc Markdown.doc"""
    newton(R::QQMPolyRing, bb::Function, d::Vector{Int}; <keyword arguments>)

Compute the multivariate polynomial corresponding to the black box function `bb` of degrees at most `d` using Newton's interpolation formula.

# Arguments
- `R::QQMPolyRing`: the multivariate polynomial ring over the rationals.
- `bb::Function`: a black box function that takes a vector of rational numbers as input and returns a rational number as output.
- `d::Vector{Int}`: the maximum degrees of the interpolating polynomial in each variable.
- `non_vanishing_poly::QQMPolyRingElem=one(R)`: a polynomial that does not vanish on the points to be evaluated. This can be used to avoid evaluating the black box function at points where it is not defined.
- `nr_thrds::Int=1`: the number of threads to use when evaluating the black box function.
- `show_progress::Bool=false`: whether to show a progress bar while collecting points.
- `desc::String="Multivariate polynomial interpolation"`: the description to show in the progress bar.
"""
function newton(
    R::QQMPolyRing,
    bb::Function,
    d::Vector{Int};
    non_vanishing_poly::QQMPolyRingElem=one(R),
    nr_thrds::Int=1,
    show_progress::Bool=false,
    desc::String="Multivariate polynomial interpolation"
)::QQMPolyRingElem
    @assert !iszero(non_vanishing_poly)
    t = length(gens(R))
    @assert length(d) == t
    x = Vector{Vector{ZZRingElem}}()
    y = Vector{QQFieldElem}()
    data_lock = ReentrantLock()
    prog = ProgressBar(total=prod(d .+ 1); desc=desc, enabled=show_progress)
    update!(prog, 0)
    function populate(cur::Vector{ZZRingElem}, dim::Int; num_threads::Int=1)
        if dim > t
            lock(data_lock) do
                push!(x, cur)
                push!(y, bb(QQ.(cur)))
                next!(prog)
            end
            return
        end
        total = Atomic(0)
        i = 1
        while total.x < d[dim] + 1
            num_threads_chunk = min(num_threads, d[dim] + 1 - total.x)
            Threads.@threads for j in 0:num_threads_chunk-1
                if !iszero(evaluate(non_vanishing_poly, gens(R)[1:length(cur)+1], [cur; ZZ(i + j)]))
                    populate([cur; ZZ(i + j)], dim + 1)
                    @atomic total.x += 1
                end
            end
            i += num_threads_chunk
        end
    end
    populate(Vector{ZZRingElem}(), 1; num_threads=nr_thrds)
    f = newton(R, x, y, d)
    finish!(prog)
    return f
end
