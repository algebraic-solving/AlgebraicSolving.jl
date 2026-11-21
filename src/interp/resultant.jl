# TODO: remove this once the next version of Nemo is released
function _flint_discriminant(f::QQMPolyRingElem, i::Int)::QQMPolyRingElem
    R = parent(f)
    res = zero(R)
    @ccall Nemo.libflint.fmpq_mpoly_discriminant(res::Ref{QQMPolyRingElem}, f::Ref{QQMPolyRingElem}, (i - 1)::Int, R::Ref{QQMPolyRing})::Int
    return res
end

# TODO: remove this once the next version of Nemo is released
function _flint_resultant(f::QQMPolyRingElem, g::QQMPolyRingElem, i::Int)::QQMPolyRingElem
    R = parent(f)
    res = zero(R)
    @ccall Nemo.libflint.fmpq_mpoly_resultant(res::Ref{QQMPolyRingElem}, f::Ref{QQMPolyRingElem}, g::Ref{QQMPolyRingElem}, (i - 1)::Int, R::Ref{QQMPolyRing})::Int
    return res
end

# Map p to the ring R by sending the k-th variable of parent(p) to the c[k]-th variable of R, or to zero if c[k] == 0.
function _flint_map_gens(R::QQMPolyRing, p::QQMPolyRingElem, c::Vector{Int})::QQMPolyRingElem
    S = parent(p)
    @assert length(c) == nvars(S) && all(k -> 0 <= k <= nvars(R), c)
    c0 = Int[k - 1 for k in c]
    z = zero(R)
    @ccall Nemo.libflint.fmpq_mpoly_compose_fmpq_mpoly_gen(z::Ref{QQMPolyRingElem}, p::Ref{QQMPolyRingElem}, c0::Ptr{Int}, S::Ref{QQMPolyRing}, R::Ref{QQMPolyRing})::Nothing
    return z
end

# The univariate polynomial f(vals[1], ..., vals[i-1], t, vals[i], ...) in R.
function _flint_specialize(R::QQPolyRing, f::QQMPolyRingElem, i::Int, vals::Vector{QQFieldElem})::QQPolyRingElem
    S = parent(f)
    n = nvars(S)
    @assert length(vals) == n - 1 && 1 <= i <= n
    g = _flint_evaluate(f, [1:i-1; i+1:n], vals)
    z = R()
    ok = @ccall Nemo.libflint.fmpq_mpoly_get_fmpq_poly(z::Ref{QQPolyRingElem}, g::Ref{QQMPolyRingElem}, (i - 1)::Int, S::Ref{QQMPolyRing})::Cint
    @assert ok != 0 "fmpq_mpoly_get_fmpq_poly failed"
    return z
end

function discriminant(
    f::QQMPolyRingElem,
    i::Int;
    nr_thrds::Int=1,
    show_progress::Bool=false,
    desc::String="Discriminant via interpolation"
)::QQMPolyRingElem
    R = parent(f)
    res = zero(R)
    R′, _ = polynomial_ring(QQ, [symbols(R)[1:i-1]; symbols(R)[i+1:end]])
    others = [1:i-1; i+1:nvars(R)]
    to_R′ = [k < i ? k : (k == i ? 0 : k - 1) for k in 1:nvars(R)]
    to_R = others
    if degree(f, i) >= 1
        Qt, _ = polynomial_ring(QQ, :t; cached=false)
        n = degree(f, i)
        sgn = isodd(div(n * (n - 1), 2)) ? -1 : 1
        bb = function (v)
            p = _flint_specialize(Qt, f, i, v)
            @assert degree(p) == n
            return sgn * Nemo.resultant(p, derivative(p)) // leading_coefficient(p)
        end
    else
        bb = v -> constant_coefficient(_flint_discriminant(_flint_evaluate(f, others, v), i))
    end
    deg_f = total_degree(f)
    d = deg_f * (deg_f - 1)
    lc = _flint_map_gens(R′, leading_coefficient(f, i), to_R′)
    res = newton(R′, bb, d, non_vanishing_poly=lc, nr_thrds=nr_thrds, show_progress=show_progress, desc=desc)
    @assert total_degree(res) <= d
    res = _flint_map_gens(R, res, to_R)
    return res
end

function resultant(
    f::QQMPolyRingElem,
    g::QQMPolyRingElem,
    i::Int;
    nr_thrds::Int=1,
    show_progress::Bool=false,
    desc::String="Resultant via interpolation"
)::QQMPolyRingElem
    @assert parent(f) == parent(g)
    R = parent(f)
    res = zero(R)
    R′, _ = polynomial_ring(QQ, [symbols(R)[1:i-1]; symbols(R)[i+1:end]])
    others = [1:i-1; i+1:nvars(R)]
    to_R′ = [k < i ? k : (k == i ? 0 : k - 1) for k in 1:nvars(R)]
    to_R = others
    if degree(f, i) >= 1 && degree(g, i) >= 1
        Qt, _ = polynomial_ring(QQ, :t; cached=false)
        bb = v -> Nemo.resultant(_flint_specialize(Qt, f, i, v), _flint_specialize(Qt, g, i, v))
    else
        bb = v -> constant_coefficient(_flint_resultant(_flint_evaluate(f, others, v), _flint_evaluate(g, others, v), i))
    end
    deg_f = total_degree(f)
    deg_g = total_degree(g)
    d = deg_f * deg_g
    lc = _flint_map_gens(R′, leading_coefficient(f, i) * leading_coefficient(g, i), to_R′)
    res = newton(R′, bb, d, non_vanishing_poly=lc, nr_thrds=nr_thrds, show_progress=show_progress, desc=desc)
    @assert total_degree(res) <= d
    res = _flint_map_gens(R, res, to_R)
    return res
end
