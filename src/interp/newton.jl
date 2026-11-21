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
    n == 0 && return z
    ax = Int[x[i].d for i in 1:n]
    ay = Fmpq[Fmpq(y[i].num, y[i].den) for i in 1:n]
    ok = GC.@preserve x y begin
        @ccall Nemo.libflint.fmpq_poly_interpolate_fmpz_fmpq_vec(z::Ref{QQPolyRingElem}, ax::Ptr{Int}, ay::Ptr{Fmpq}, n::Int)::Cint
    end
    @assert ok != 0 "interpolation nodes must be pairwise distinct"
    @assert evaluate(z, x[n]) == y[n]
    return z
end

@doc Markdown.doc"""
    newton(R::QQMPolyRing, x::Vector{Vector{ZZRingElem}}, y::Vector{QQFieldElem}, d::Vector{Int})

Compute the multivariate polynomial corresponding to the given points `x` and values `y` of degrees at most `d` using Newton's interpolation formula.

**Note**: The points in `x` should form a grid of size `(d[1]+1) x (d[2]+1) x ... x (d[n]+1)`, where `n` is the number of variables.
The points should be ordered such that the first variable changes the slowest and the last variable changes the fastest.
For example, for `d = [1, 1]`, the points should be ordered as `[[0, 0], [0, 1], [1, 0], [1, 1]]`.

# Arguments
- `R::QQMPolyRing`: the multivariate polynomial ring over the rationals.
- `x::Vector{Vector{ZZRingElem}}`: the points at which the polynomial is evaluated. Each point is a vector of rational numbers represented as `ZZRingElem`.
- `y::Vector{QQFieldElem}`: the values of the polynomial at the points `x`.
- `d::Vector{Int}`: the maximum degrees of the interpolating polynomial in each variable.

# Example
```jldoctest
julia> using AlgebraicSolving

julia> using AlgebraicSolving.Interpolation

julia> R, _ = polynomial_ring(QQ, ["a", "b"]);

julia> x = Vector{ZZRingElem}[[0, 0], [0, 1], [1, 0], [1, 1]];

julia> y = QQFieldElem[0, 3, 2, 5];

julia> d = [1, 1];

julia> f = newton(R, x, y, d)
2*a + 3*b
```

**Note**: This is an internal function.
"""
function newton(
    R::QQMPolyRing,
    x::Vector{Vector{ZZRingElem}},
    y::Vector{QQFieldElem},
    d::Vector{Int}
)::QQMPolyRingElem
    @assert length(x) == length(y) == prod(d .+ 1)
    if isempty(d)
        return R(y[1])
    end
    n = length(x)
    step = prod(d[2:end] .+ 1)
    x_index = length(gens(R)) - length(d) + 1
    z = gens(R)[x_index]
    x_trunc = [x[i][x_index] for i in 1:step:n]
    for i in 1:n
        @assert x[i][x_index] == x_trunc[(i-1)÷step+1]
    end
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

# Example
```jldoctest
julia> using AlgebraicSolving

julia> using AlgebraicSolving.Interpolation

julia> R, _ = polynomial_ring(QQ, ["a", "b"]);

julia> f = t -> (t[1] + t[2])^2 + (t[1] - t[2])^2;

julia> d = [2, 2];

julia> newton(R, f, d)
2*a^2 + 2*b^2
```
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
    perm = sortperm(x)
    x = x[perm]
    y = y[perm]
    f = newton(R, x, y, d)
    finish!(prog)
    return f
end

@doc Markdown.doc"""
    _simplex_indices(n::Int, d::Int)

Return all multi-indices `a` in $\mathbb{N}^n$ with `sum(a) <= d`, sorted first by total degree and then lexicographically.
In particular, for every `i <= d`, the multi-indices with `sum(a) <= i` are exactly the first `binomial(n + i, n)` entries.

**Note**: This is an internal function.
"""
function _simplex_indices(n::Int, d::Int)::Vector{Vector{Int}}
    @assert n >= 0 && d >= 0
    res = Vector{Vector{Int}}()
    if n == 0
        push!(res, Int[])
        return res
    end
    cur = zeros(Int, n)
    # all a in N^n with sum(a) == s, in lexicographically decreasing order of a[1], a[2], ...
    function rec(j::Int, s::Int)
        if j == n
            cur[n] = s
            push!(res, copy(cur))
            return
        end
        for k in s:-1:0
            cur[j] = k
            rec(j + 1, s - k)
        end
    end
    for s in 0:d
        rec(1, s)
    end
    return res
end

@doc Markdown.doc"""
    _flint_change_basis!(p::QQPolyRingElem, roots::Vector{ZZRingElem}, to_newton::Bool)

Convert the coefficients of `p` in place from the monomial basis to the Newton basis
$1, (t - r_0), (t - r_0)(t - r_1), \ldots$ for the given `roots` (if `to_newton`), or back.

**Note**: This is an internal function.
"""
function _flint_change_basis!(
    p::QQPolyRingElem,
    roots::Vector{ZZRingElem},
    to_newton::Bool
)::QQPolyRingElem
    len = length(p)
    len <= 1 && return p
    @assert length(roots) >= len - 1
    # contiguous fmpz handles, see _flint_lagrange
    ar = Int[roots[i].d for i in 1:len-1]
    GC.@preserve roots p begin
        if to_newton
            @ccall Nemo.libflint._fmpz_poly_monomial_to_newton(p.coeffs::Ptr{ZZRingElem}, ar::Ptr{Int}, len::Int)::Nothing
        else
            @ccall Nemo.libflint._fmpz_poly_newton_to_monomial(p.coeffs::Ptr{ZZRingElem}, ar::Ptr{Int}, len::Int)::Nothing
        end
        # both changes of basis are unitriangular, so the length is unchanged;
        # only the gcd with the denominator may have changed
        @ccall Nemo.libflint.fmpq_poly_canonicalise(p::Ref{QQPolyRingElem})::Nothing
    end
    return p
end

@doc Markdown.doc"""
    _newton_simplex(R::QQMPolyRing, var_idx::Vector{Int}, nodes::Vector{Vector{ZZRingElem}}, y::Vector{QQFieldElem}, d::Int)

Compute the polynomial of total degree at most `d` in the variables `gens(R)[var_idx]` taking the values `y` on the simplex lattice built from the per-axis `nodes`.

The value `y[k]` must be the value at the point `[nodes[j][a[j] + 1] for j in eachindex(var_idx)]`, where `a = _simplex_indices(length(var_idx), d)[k]`.
Only the first `d + 1` nodes of each axis are used, and these must be pairwise distinct.

**Note**: This is an internal function.
"""
function _newton_simplex(
    R::QQMPolyRing,
    var_idx::Vector{Int},
    nodes::Vector{Vector{ZZRingElem}},
    y::Vector{QQFieldElem},
    d::Int
)::QQMPolyRingElem
    n = length(var_idx)
    @assert length(nodes) == n
    @assert all(j -> length(nodes[j]) >= d + 1 && allunique(nodes[j][1:d+1]), 1:n)
    idx = _simplex_indices(n, d)
    @assert length(y) == length(idx)
    if n == 0
        return R(y[1])
    end
    Qt, _ = polynomial_ring(QQ, :t; cached=false)
    if n == 1
        p = _flint_lagrange(Qt, nodes[1][1:d+1], y)
        c = [coeff(p, k) for k in 0:d]
    else
        pos = Dict{Vector{Int},Int}(a => k for (k, a) in enumerate(idx))
        c = copy(y)
        line = Vector{Int}(undef, d + 1)
        # Collect the positions of the lines a, a + e_j, ..., a + L e_j with a[j] == 0.
        function _line!(a::Vector{Int}, j::Int)::Int
            L = d - sum(a)
            b = copy(a)
            for k in 0:L
                b[j] = k
                line[k+1] = pos[b]
            end
            return L
        end
        # 1. Newton coefficients along each axis. The first L + 1 Newton coefficients of
        #    a univariate function only depend on its values at the first L + 1 nodes, so
        #    they are the Newton coefficients of the degree-L interpolant of the line.
        for j in 1:n
            x = nodes[j][1:d+1]
            for a in idx
                a[j] == 0 || continue
                L = _line!(a, j)
                L == 0 && continue
                p = _flint_lagrange(Qt, x[1:L+1], c[line[1:L+1]])
                _flint_change_basis!(p, x, true)
                for k in 0:L
                    c[line[k+1]] = coeff(p, k)
                end
            end
        end
        # 2. Change of basis from Newton to monomials along each axis.
        for j in 1:n
            x = nodes[j][1:d+1]
            for a in idx
                a[j] == 0 || continue
                L = _line!(a, j)
                L == 0 && continue
                p = Qt(c[line[1:L+1]])
                _flint_change_basis!(p, x, false)
                for k in 0:L
                    c[line[k+1]] = coeff(p, k)
                end
            end
        end
    end
    # 3. Assemble the polynomial.
    C = MPolyBuildCtx(R)
    exp = zeros(Int, length(gens(R)))
    for (k, a) in enumerate(idx)
        iszero(c[k]) && continue
        for j in 1:n
            exp[var_idx[j]] = a[j]
        end
        push_term!(C, c[k], exp)
    end
    return finish(C)
end

@doc Markdown.doc"""
    _simplex_point(nodes::Vector{Vector{ZZRingElem}}, a::Vector{Int})

Return the point of the simplex lattice with multi-index `a`.

**Note**: This is an internal function.
"""
_simplex_point(nodes::Vector{Vector{ZZRingElem}}, a::Vector{Int})::Vector{ZZRingElem} =
    [nodes[j][a[j]+1] for j in eachindex(a)]

@doc Markdown.doc"""
    _replace_node!(nodes::Vector{Vector{ZZRingElem}}, a::Vector{Int})

Replace one node used by the lattice point with multi-index `a` by a fresh integer larger than all current nodes, and return the pair `(j, a[j])` identifying the replaced node.
The axis `j` is the one where `a[j]` is the largest, so that the lattice points sharing the replaced node are as few as possible.

**Note**: This is an internal function.
"""
function _replace_node!(
    nodes::Vector{Vector{ZZRingElem}},
    a::Vector{Int}
)::Tuple{Int,Int}
    @assert !isempty(a)
    j = argmax([(a[j], j) for j in eachindex(a)])
    nodes[j][a[j]+1] = maximum(maximum, nodes) + 1
    return j, a[j]
end

# FLINT wrapper replacing the generic AbstractAlgebra method evaluate(f, vars, vals).
function _flint_evaluate(f::QQMPolyRingElem, idx::Vector{Int}, vals::Vector{QQFieldElem})::QQMPolyRingElem
    @assert length(idx) == length(vals)
    R = parent(f)
    res = f
    for (j, v) in zip(idx, vals)
        @assert 1 <= j <= nvars(R)
        z = zero(R)
        ok = @ccall Nemo.libflint.fmpq_mpoly_evaluate_one_fmpq(z::Ref{QQMPolyRingElem}, res::Ref{QQMPolyRingElem}, (j - 1)::Int, v::Ref{QQFieldElem}, R::Ref{QQMPolyRing})::Cint
        # FLINT returns 0 only if the result cannot be computed (exponent/size overflow)
        @assert ok != 0 "fmpq_mpoly_evaluate_one_fmpq failed"
        res = z
    end
    return res
end

@doc Markdown.doc"""
    _choose_nodes(p::QQMPolyRingElem, d::Int)

Return per-axis nodes such that the nonzero polynomial `p` does not vanish at any point of the simplex lattice of total degree `d` built from them.
If `p` does not vanish on the default lattice, the nodes are `1, ..., d + 1` on every axis.

**Note**: This is an internal function.
"""
function _choose_nodes(p::QQMPolyRingElem, d::Int)::Vector{Vector{ZZRingElem}}
    @assert !iszero(p) && d >= 0
    R = parent(p)
    n = nvars(R)
    is_constant(p) && return [ZZ.(1:d+1) for _ in 1:n]
    nodes = [ZZRingElem[] for _ in 1:n]
    # partial evaluations of p at the prefixes (a_1, ..., a_{j-1}) built so far
    prefixes = Dict{Vector{Int},QQMPolyRingElem}(Int[] => p)
    for j in 1:n
        next_prefixes = Dict{Vector{Int},QQMPolyRingElem}()
        c = 0  # candidates only increase, so the nodes of an axis are distinct
        for k in 0:d
            relevant = [(a, q) for (a, q) in prefixes if sum(a) + k <= d]
            while true
                c += 1
                vals = [_flint_evaluate(q, [j], [QQ(c)]) for (_, q) in relevant]
                if all(!iszero, vals)
                    push!(nodes[j], ZZ(c))
                    for ((a, _), q) in zip(relevant, vals)
                        next_prefixes[[a; k]] = q
                    end
                    break
                end
            end
        end
        prefixes = next_prefixes
    end
    return nodes
end

@doc Markdown.doc"""
    newton(R::QQMPolyRing, bb::Function, d::Int; <keyword arguments>)

Compute the multivariate polynomial of total degree at most `d` corresponding to the black box function `bb` using Newton's interpolation formula on a simplex lattice.

# Arguments
- `R::QQMPolyRing`: the multivariate polynomial ring over the rationals.
- `bb::Function`: a black box function that takes a vector of rational numbers as input and returns a rational number as output.
- `d::Int`: an upper bound on the total degree of the interpolating polynomial.
- `non_vanishing_poly::QQMPolyRingElem=one(R)`: a polynomial that does not vanish on the points to be evaluated. This can be used to avoid evaluating the black box function at points where it is not defined.
- `nr_thrds::Int=1`: the number of threads to use when evaluating the black box function.
- `show_progress::Bool=false`: whether to show a progress bar while collecting points.
- `desc::String="Multivariate polynomial interpolation"`: the description to show in the progress bar.

# Example
```jldoctest
julia> using AlgebraicSolving

julia> using AlgebraicSolving.Interpolation

julia> R, _ = polynomial_ring(QQ, ["a", "b"]);

julia> f = t -> (t[1] + t[2])^2 + (t[1] - t[2])^2;

julia> newton(R, f, 2)
2*a^2 + 2*b^2
```
"""
function newton(
    R::QQMPolyRing,
    bb::Function,
    d::Int;
    non_vanishing_poly::QQMPolyRingElem=one(R),
    nr_thrds::Int=1,
    show_progress::Bool=false,
    desc::String="Multivariate polynomial interpolation"
)::QQMPolyRingElem
    @assert d >= 0
    @assert !iszero(non_vanishing_poly)
    @assert parent(non_vanishing_poly) == R
    n = length(gens(R))
    nodes = _choose_nodes(non_vanishing_poly, d)
    idx = _simplex_indices(n, d)
    N = length(idx)
    y = Vector{QQFieldElem}(undef, N)
    prog = ProgressBar(total=N; desc=desc, enabled=show_progress)
    update!(prog, 0)
    prog_lock = ReentrantLock()
    function eval_range(r)
        for k in r
            y[k] = bb(QQ.(_simplex_point(nodes, idx[k])))
            lock(() -> next!(prog), prog_lock)
        end
    end
    nt = max(1, min(nr_thrds, N))
    if nt == 1
        eval_range(1:N)
    else
        chunks = [k:nt:N for k in 1:nt]
        Threads.@threads for r in chunks
            eval_range(r)
        end
    end
    f = _newton_simplex(R, collect(1:n), nodes, y, d)
    finish!(prog)
    return f
end
