function _map_exponent_vectors(f, p::MPolyRingElem, R::MPolyRing)::MPolyRingElem
    cvzip = zip(coefficients(p), exponent_vectors(p))
    M = MPolyBuildCtx(R)
    for (c, v) in cvzip
        push_term!(M, coefficient_ring(R)(c), f(v))
    end
    finish(M)
end

function _change_ring(p::MPolyRingElem, R::MPolyRing)::MPolyRingElem
    n = number_of_variables(parent(p))
    n′ = number_of_variables(R)
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

function _change_ring(v::Vector{<:MPolyRingElem}, R::MPolyRing)::Vector{<:MPolyRingElem}
    map(p -> _change_ring(p, R), v)
end

function _change_ring(A::MatSpaceElem{<:MPolyRingElem}, R::MPolyRing)::MatSpaceElem{<:MPolyRingElem}
    map(f -> _change_ring(f, R), A)
end

function _saturate(
    f::Vector{QQMPolyRingElem},
    g::QQMPolyRingElem;
    worker_pool::AbstractWorkerPool=default_worker_pool()
)::Vector{QQMPolyRingElem}
    @assert !isempty(f)
    R = parent(f[1])
    n = number_of_variables(R)
    R′, vars = polynomial_ring(QQ, [:sat; symbols(R)])
    f′ = _change_ring(f, R′)
    g′ = _change_ring(g, R′)
    gb = AlgebraicSolving.groebner_basis(
        AlgebraicSolving.Ideal([f′; vars[1] * g′ - 1]);
        eliminate=1,
        worker_pool=worker_pool
    )
    gb = _change_ring(gb, R)
    gb
end

function _radical(
    f::Vector{QQMPolyRingElem};
    worker_pool::AbstractWorkerPool=default_worker_pool()
)::Vector{QQMPolyRingElem}
    @assert !isempty(f)
    R = parent(f[1])
    n = number_of_variables(R)
    s = symbols(R)
    rads = Vector{QQMPolyRingElem}()
    for _ in 1:n
        R′, _ = polynomial_ring(QQ, s, internal_ordering=internal_ordering(R))
        disc = AlgebraicSolving.groebner_basis(
            AlgebraicSolving.Ideal(_change_ring(f, R′));
            eliminate=n - 1,
            worker_pool=worker_pool
        )
        r = R′(1)
        for (p, _) in factor(prod(disc))
            r *= p
        end
        push!(rads, _change_ring(r, R))
        s = [s[2:end]; s[1]]
    end
    gb = AlgebraicSolving.groebner_basis(
        AlgebraicSolving.Ideal([f; rads]);
        worker_pool=worker_pool
    )
    gb = _change_ring(gb, R)
    gb
end

function groebner_basis(
    I::ParametricIdeal{K},
    vals::Vector{QQFieldElem};
    worker_pool::AbstractWorkerPool=default_worker_pool()
)::Vector{QQMPolyRingElem} where {K<:FracFieldElem}
    gb_spec_hash = hash((I.gens, vals))
    if haskey(I.gb_spec, gb_spec_hash)
        return I.gb_spec[gb_spec_hash]
    end
    n = I.num_vars
    elim = number_of_variables(parent(I.gens[1])) - n
    R = parent(I.gens[1])
    R₀, _ = polynomial_ring(QQ, symbols(R)[1:end], internal_ordering=:degrevlex)
    R₁, _ = polynomial_ring(QQ, I.symbols, internal_ordering=:degrevlex)
    fₓ = _change_ring(map(p -> map_coefficients(k -> evaluate(k, vals), p), I.gens), R₀)
    gb = AlgebraicSolving.groebner_basis(
        AlgebraicSolving.Ideal(fₓ);
        eliminate=elim,
        worker_pool=worker_pool
    )
    gb = _change_ring(gb, R₁)
    for sat in I.sats
        gb = _saturate(gb, _change_ring(map_coefficients(k -> evaluate(k, vals), sat), R₁); worker_pool=worker_pool)
    end
    if I.radicalize
        gb = _radical(gb; worker_pool=worker_pool)
    end
    gb = map(p -> p / leading_coefficient(p), gb)
    lock(I.gb_spec_lock) do
        I.gb_spec[gb_spec_hash] = gb
    end
    gb
end

function _is_zero_dimensional(gb::Vector{T})::Bool where {T<:MPolyRingElem}
    if isempty(gb)
        return false
    end
    R = parent(gb[1])
    n = number_of_variables(R)
    has_power = falses(n)
    for f in gb
        indices = var_indices(Nemo.leading_monomial(f))
        if isempty(indices)
            return true
        end
        if length(indices) == 1
            has_power[indices[1]] = true
        end
    end
    all(has_power)
end

function _monomial_basis(gb::Vector{T})::Vector{T} where {T<:MPolyRingElem}
    if !_is_zero_dimensional(gb)
        error("The ideal is not zero-dimensional, so a finite monomial basis does not exist.")
    end
    R = parent(gb[1])
    n = number_of_variables(R)
    lms = [Nemo.leading_monomial(f) for f in gb]
    basis = typeof(gb[1])[]
    queue = [one(R)]
    while !isempty(queue)
        u = popfirst!(queue)
        if all(!divides(u, lm)[1] for lm in lms) && !(u in basis)
            push!(basis, u)
            for i in n:-1:1
                push!(queue, u * gens(R)[i])
            end
        end
    end
    basis
end

function _expected_degree(
    I::ParametricIdeal{K};
    worker_pool::AbstractWorkerPool=default_worker_pool()
)::Int where {K<:FracFieldElem}
    gb = groebner_basis(I, QQ.(I.shift); worker_pool=worker_pool)
    b = _monomial_basis(gb)
    length(b)
end

@doc Markdown.doc"""
    groebner_basis(I::ParametricIdeal{K}; <keyword arguments>) -> Vector{Generic.MPoly{FracFieldElem{QQMPolyRingElem}}}

Computes a parametric Gröbner basis of the parametric ideal `I` w.r.t. the degree reverse lexicographic ordering. The Gröbner basis is computed by evaluating the ideal at a number of parameter values and interpolating the coefficients of the resulting Gröbner bases.

# Arguments
- `retry::Int=10`: the maximum number of consecutive failures allowed when computing the Gröbner basis or interpolating the coefficients.
- `nr_thrds::Int=1`: the number of threads to use for parallel computations.
- `worker_pool::AbstractWorkerPool=default_worker_pool()`: the worker pool to use for parallel computations.
- `show_progress::Bool=false`: whether to show a progress bar while computing the Gröbner basis.
"""
function groebner_basis(
    I::ParametricIdeal{K};
    retry::Int=10,
    nr_thrds::Int=1,
    worker_pool::AbstractWorkerPool=default_worker_pool(),
    show_progress::Bool=false
)::Vector{Generic.MPoly{FracFieldElem{QQMPolyRingElem}}} where {K<:FracFieldElem}
    input_hash = hash(I.gens)
    if haskey(I.gb, input_hash)
        return I.gb[input_hash]
    end
    R = parent(I.gens[1])
    R₂ = base_ring(I.base_frac_field)
    R₄, _ = polynomial_ring(I.base_frac_field, I.symbols, internal_ordering=internal_ordering(R))
    gb_shift_mons = nothing
    function bb(v)
        gb = groebner_basis(I, v; worker_pool=worker_pool)
        if !isnothing(gb_shift_mons)
            gb_mons = [collect(monomials(p)) for p in gb]
            if gb_mons != gb_shift_mons
                error("Error evaluating black box function: the Gröbner basis has unexpected monomial structure at parameter values $(v).")
            end
        end
        gb
    end
    gb_shift = bb(QQ.(I.shift))
    gb_shift_mons = [collect(monomials(p)) for p in gb_shift]
    pgb = fill(R₄(0), length(gb_shift))
    len = length(gb_shift)
    len_mons = length.(gb_shift)
    for i in 1:len
        for j in 1:len_mons[i]
            denom = R₂(1)
            c = Interpolation.cuyt_lee(
                R₂, v -> evaluate(denom, v) * coeff(bb(v)[i], j);
                initial_shift=I.shift,
                retry=retry,
                nr_thrds=nr_thrds,
                show_progress=show_progress,
                desc="Parametric Gröbner basis, size=($len,$(maximum(len_mons))), index=($i,$j)"
            )
            pgb[i] += R₄(c) / denom * _change_ring(Nemo.monomial(gb_shift[i], j), R₄)
        end
    end
    lock(I.gb_lock) do
        I.gb[input_hash] = pgb
    end
    pgb
end
