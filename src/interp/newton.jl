struct Fmpq
    num::Int
    den::Int
end

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

function _from_univariate(
    x::QQMPolyRingElem,
    f::QQPolyRingElem
)::QQMPolyRingElem
    R = parent(x)
    j = findfirst(u -> u == x, gens(R))
    t = length(gens(R))
    C = MPolyBuildCtx(R)
    for i in 0:length(f)-1
        c = coeff(f, i)
        if c != 0
            exp = fill(0, t)
            exp[j] = i
            push_term!(C, c, exp)
        end
    end
    finish(C)
end

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
        return _from_univariate(z, f)
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
