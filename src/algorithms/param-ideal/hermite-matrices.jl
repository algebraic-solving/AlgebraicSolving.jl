@doc Markdown.doc"""
    hermite_matrix(gb::Vector{T}, b::Vector{T}, g::T) -> MatElem

Computes the Hermite matrix associated with the polynomial `g` with respect to the Groebner basis `gb` and the monomial basis `b`.
"""
function hermite_matrix(gb::Vector{T}, b::Vector{T}, g::T)::MatElem where {T<:MPolyRingElem}
    δ = length(b)
    H = zero_matrix(coefficient_ring(gb[1]), δ, δ)
    M_g = multiplication_matrix(gb, b, g)
    M = [multiplication_matrix(gb, b, b[i]) for i in 1:δ]
    for i in 1:δ
        for j in i:δ
            H[i, j] = tr(M[i] * M[j] * M_g)
        end
    end
    for i in 2:δ
        for j in 1:(i-1)
            H[i, j] = H[j, i]
        end
    end
    H
end

@doc Markdown.doc"""
    hermite_matrix(I::ParametricIdeal{K}, g::MPoly{K}; <keyword arguments>) -> MatSpaceElem{K}

Computes the parametric Hermite matrix associated with the polynomial `g` with respect to the parametric ideal `I`. The Hermite matrix is computed by evaluating the ideal at a number of parameter values and interpolating the coefficients of the resulting Hermite matrices.

# Arguments
- `retry::Int=10`: the maximum number of consecutive failures allowed when computing the Hermite matrix or interpolating the coefficients.
- `nr_thrds::Int=1`: the number of threads to use for parallel computations.
- `worker_pool::AbstractWorkerPool=default_worker_pool()`: the worker pool to use for parallel computations.
- `show_progress::Bool=false`: whether to show a progress bar while computing the Hermite matrix.
"""
function hermite_matrix(
    I::ParametricIdeal{K},
    g::MPoly{K};
    retry::Int=10,
    nr_thrds::Int=1,
    worker_pool::AbstractWorkerPool=default_worker_pool(),
    show_progress::Bool=false
)::MatSpaceElem{K} where {K<:FracFieldElem}
    input_hash = hash((I.gens, g))
    if haskey(I.hm, input_hash)
        return I.hm[input_hash]
    end
    if show_progress
        @info "Computing parametric Hermite matrix"
    end
    if I.prefer_gb
        gb = groebner_basis(I; retry=retry, nr_thrds=nr_thrds, worker_pool=worker_pool, show_progress=show_progress)
        b = _monomial_basis(gb)
        g′ = _change_ring(g, parent(gb[1]))
        H = hermite_matrix(gb, b, g′)
        H = map(_canonicalize, H)
        lock(I.hm_lock) do
            I.hm[input_hash] = H
        end
        return H
    end
    δ = _expected_degree(I; worker_pool=worker_pool)
    R₁, _ = polynomial_ring(QQ, I.symbols, internal_ordering=:degrevlex)
    R₂ = base_ring(I.base_frac_field)
    function bb(v)
        hm_spec_hash = hash((I.gens, g, v))
        if haskey(I.hm_spec, hm_spec_hash)
            return I.hm_spec[hm_spec_hash]
        end
        gb = groebner_basis(I, v; worker_pool=worker_pool)
        b = _monomial_basis(gb)
        if length(b) != δ
            error("Error evaluating black box function: the ideal has unexpected degree at parameter values $(v).")
        end
        gₓ = _change_ring(map_coefficients(k -> evaluate(k, v), g), R₁)
        Hₓ = hermite_matrix(gb, b, gₓ)
        lock(I.hm_spec_lock) do
            I.hm_spec[hm_spec_hash] = Hₓ
        end
        Hₓ
    end
    H = zero_matrix(fraction_field(R₂), δ, δ)
    for i in 1:δ
        for j in i:δ
            H[i, j] = Interpolation.cuyt_lee(
                R₂, v -> bb(v)[i, j];
                initial_shift=I.shift,
                retry=retry,
                nr_thrds=nr_thrds,
                show_progress=show_progress,
                desc="Parametric Hermite matrix, size=($δ,$δ), index=($i,$j)"
            )
        end
    end
    for i in 2:δ
        for j in 1:(i-1)
            H[i, j] = H[j, i]
        end
    end
    lock(I.hm_lock) do
        I.hm[input_hash] = H
    end
    H
end

# To compute Her(g^α), computing as above is a bit wasteful
# We proceed by first computing Her(1) and then multiplying by Mul(g)^α as needed
@doc Markdown.doc"""
    hermite_matrix(g::Vector{MPoly{K}}, α::Vector{Int}, I::ParametricIdeal{K}; <keyword arguments>) -> MatSpaceElem{K}
Computes the parametric Hermite matrix associated with the polynomials `g` raised to the powers specified in `α` with respect to the parametric ideal `I`. The Hermite matrix is computed by evaluating the ideal at a number of parameter values and interpolating the coefficients of the resulting Hermite matrices.

# Arguments
- `nr_thrds::Int=1`: the number of threads to use for parallel computations.
- `worker_pool::AbstractWorkerPool=default_worker_pool()`: the worker pool to use for parallel computations.
- `show_progress::Bool=false`: whether to show a progress bar while computing the Hermite matrix.
"""
function hermite_matrix(
    g::Vector{MPoly{K}},
    α::Vector{Int},
    I::ParametricIdeal{K};
    nr_thrds::Int=1,
    worker_pool::AbstractWorkerPool=default_worker_pool(),
    show_progress::Bool=false
)::MatSpaceElem{K} where {K<:FracFieldElem}
    R = parent(I.gens[1])
    H = hermite_matrix(I, one(R); nr_thrds=nr_thrds, worker_pool=worker_pool, show_progress=show_progress)
    for i in 1:length(g)
        if α[i] == 0
            continue
        end
        M = multiplication_matrix(I, g[i]; nr_thrds=nr_thrds, worker_pool=worker_pool, show_progress=show_progress)
        if α[i] == 1
            H = H * M
        else
            @assert α[i] == 2
            H = H * M * M
        end
    end
    H = map(_canonicalize, H)
    H
end

@doc Markdown.doc"""
    hermite_matrix(g::Vector{MPoly{K}}, α::Vector{Int}, I::ParametricIdeal{K}, vals::Vector{QQFieldElem}; <keyword arguments>) -> QQMatrix
Computes the Hermite matrix associated with the polynomials `g` raised to the powers specified in `α` with respect to the parametric ideal `I` evaluated at the parameter values specified in `vals`. The Hermite matrix is computed by evaluating the ideal at the specified parameter values and interpolating the coefficients of the resulting Hermite matrices.

# Arguments
- `nr_thrds::Int=1`: the number of threads to use for parallel computations.
- `worker_pool::AbstractWorkerPool=default_worker_pool()`: the worker pool to use for parallel computations.
- `show_progress::Bool=false`: whether to show a progress bar while computing the Hermite matrix.
"""
function hermite_matrix(
    g::Vector{MPoly{K}},
    α::Vector{Int},
    I::ParametricIdeal{K},
    vals::Vector{QQFieldElem};
    nr_thrds::Int=1,
    worker_pool::AbstractWorkerPool=default_worker_pool(),
    show_progress::Bool=false
)::QQMatrix where {K<:FracFieldElem}
    R = parent(I.gens[1])
    H = map(a -> evaluate(numerator(a), vals) // evaluate(denominator(a), vals), hermite_matrix(I, one(R); nr_thrds=nr_thrds, worker_pool=worker_pool, show_progress=show_progress))
    for i in 1:length(g)
        if α[i] == 0
            continue
        end
        M = map(a -> evaluate(numerator(a), vals) // evaluate(denominator(a), vals), multiplication_matrix(I, g[i]; nr_thrds=nr_thrds, worker_pool=worker_pool, show_progress=show_progress))
        if α[i] == 1
            H = H * M
        else
            @assert α[i] == 2
            H = H * M * M
        end
    end
    H
end
