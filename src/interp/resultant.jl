function _map_exponent_vectors(f, p::MPolyRingElem, R::MPolyRing)::MPolyRingElem
    cvzip = zip(coefficients(p), exponent_vectors(p))
    M = MPolyBuildCtx(R)
    for (c, v) in cvzip
        push_term!(M, coefficient_ring(R)(c), f(v))
    end
    finish(M)
end

function _change_ring(p::MPolyRingElem, R::MPolyRing)::MPolyRingElem
    n = length(gens(parent(p)))
    n′ = length(gens(R))
    if n == 0 || n′ == 0
        return R(constant_coefficient(p))
    end
    indices = Vector{Int}()
    for s in symbols(R)
        i = findfirst(==(s), symbols(parent(p)))
        if isnothing(i)
            push!(indices, 0)
        else
            push!(indices, i)
        end
    end
    return _map_exponent_vectors(v -> [i == 0 ? 0 : v[i] for i in indices], p, R)
end

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

function discriminant(
    f::QQMPolyRingElem,
    i::Int;
    nr_thrds::Int=1,
    show_progress::Bool=false,
    desc::String="Discriminant via interpolation"
)::QQMPolyRingElem
    R = parent(f)
    res = zero(R)
    vars = [gens(R)[1:i-1]; gens(R)[i+1:end]]
    function bb(v)
        f′ = evaluate(f, vars, v)
        constant_coefficient(_flint_discriminant(f′, i))
    end
    R′, _ = polynomial_ring(QQ, [symbols(R)[1:i-1]; symbols(R)[i+1:end]])
    deg_f = total_degree(f)
    # not tight, but should be good enough in practice
    d = fill(deg_f * deg_f, length(gens(R′)))
    lc = leading_coefficient(f, i)
    res = newton(R′, bb, d, non_vanishing_poly=lc, nr_thrds=nr_thrds, show_progress=show_progress, desc=desc)
    @assert total_degree(res) <= deg_f * deg_f - deg_f
    res = _change_ring(res, R)
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
    vars = [gens(R)[1:i-1]; gens(R)[i+1:end]]
    function bb(v)
        f′ = evaluate(f, vars, v)
        g′ = evaluate(g, vars, v)
        constant_coefficient(_flint_resultant(f′, g′, i))
    end
    R′, _ = polynomial_ring(QQ, [symbols(R)[1:i-1]; symbols(R)[i+1:end]])
    deg_f = total_degree(f)
    deg_g = total_degree(g)
    # not tight, but should be good enough in practice
    d = fill(deg_f * deg_g + 1, length(gens(R′)))
    lc = leading_coefficient(f, i) * leading_coefficient(g, i)
    res = newton(R′, bb, d, non_vanishing_poly=lc, nr_thrds=nr_thrds, show_progress=show_progress, desc=desc)
    @assert total_degree(res) <= deg_f * deg_g
    res = _change_ring(res, R)
    return res
end
