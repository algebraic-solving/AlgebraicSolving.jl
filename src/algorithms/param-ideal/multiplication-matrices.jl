
@doc Markdown.doc"""
    multiplication_matrix(gb::Vector{T}, b::Vector{T}, g::T) -> MatElem

Computes the multiplication matrix associated with the polynomial `g` with respect to the Groebner basis `gb` and the monomial basis `b`.
"""
function multiplication_matrix(gb::Vector{T}, b::Vector{T}, g::T)::MatElem where {T<:MPolyRingElem}
    δ = length(b)
    M = zero_matrix(coefficient_ring(gb[1]), δ, δ)
    for i in 1:δ
        prod = divrem(b[i] * g, gb)[2]
        for j in 1:δ
            M[j, i] = coeff(prod, b[j])
        end
    end
    M
end

# For the multiplication matrices, we implement canonicalize functions
# to make computations more efficient

function _canonicalize!(p::QQMPolyRingElem)::QQMPolyRingElem
    combine_like_terms!(sort_terms!(p))
end

function _canonicalize(p::FracFieldElem{QQMPolyRingElem})::FracFieldElem{QQMPolyRingElem}
    _canonicalize!(numerator(p)) // _canonicalize!(denominator(p))
end

@doc Markdown.doc"""
    multiplication_matrix(I::ParametricIdeal{K}, g::MPoly{K}; <keyword arguments>) -> MatSpaceElem{K}

Computes the parametric multiplication matrix associated with the polynomial `g` with respect to the parametric ideal `I`. The multiplication matrix is computed by evaluating the ideal at a number of parameter values and interpolating the coefficients of the resulting multiplication matrices.

# Arguments
- `retry::Int=10`: the maximum number of consecutive failures allowed when computing the multiplication matrix or interpolating the coefficients.
- `nr_thrds::Int=1`: the number of threads to use for parallel computations.
- `worker_pool::AbstractWorkerPool=default_worker_pool()`: the worker pool to use for parallel computations.
- `show_progress::Bool=false`: whether to show a progress bar while computing the multiplication matrix.
"""
function multiplication_matrix(
    I::ParametricIdeal{K},
    g::MPoly{K};
    retry::Int=10,
    nr_thrds::Int=1,
    worker_pool::AbstractWorkerPool=default_worker_pool(),
    show_progress::Bool=false
)::MatSpaceElem{K} where {K<:FracFieldElem}
    input_hash = hash((I.gens, g))
    if haskey(I.mm, input_hash)
        return I.mm[input_hash]
    end
    if show_progress
        @info "Computing parametric multiplication matrix"
    end
    if I.prefer_gb
        gb = groebner_basis(I; retry=retry, nr_thrds=nr_thrds, worker_pool=worker_pool, show_progress=show_progress)
        b = _monomial_basis(gb)
        g′ = _change_ring(g, parent(gb[1]))
        M = multiplication_matrix(gb, b, g′)
        M = map(_canonicalize, M)
        lock(I.mm_lock) do
            I.mm[input_hash] = M
        end
        return M
    end
    δ = _expected_degree(I; worker_pool=worker_pool)
    R₁, _ = polynomial_ring(QQ, I.symbols, internal_ordering=:degrevlex)
    R₂ = base_ring(I.base_frac_field)
    function bb(v)
        mm_spec_hash = hash((I.gens, g, v))
        if haskey(I.mm_spec, mm_spec_hash)
            return I.mm_spec[mm_spec_hash]
        end
        gb = groebner_basis(I, v; worker_pool=worker_pool)
        b = _monomial_basis(gb)
        if length(b) != δ
            error("Error evaluating black box function: the ideal has unexpected degree at parameter values $(v).")
        end
        gₓ = _change_ring(map_coefficients(k -> evaluate(k, v), g), R₁)
        Mₓ = multiplication_matrix(gb, b, gₓ)
        lock(I.mm_spec_lock) do
            I.mm_spec[mm_spec_hash] = Mₓ
        end
        Mₓ
    end
    M = zero_matrix(fraction_field(R₂), δ, δ)
    for i in 1:δ
        for j in 1:δ
            M[i, j] = Interpolation.cuyt_lee(
                R₂, v -> bb(v)[i, j];
                initial_shift=I.shift,
                retry=retry,
                nr_thrds=nr_thrds,
                show_progress=show_progress,
                desc="Parametric multiplication matrix, size=($δ,$δ), index=($i,$j)"
            )
        end
    end
    lock(I.mm_lock) do
        I.mm[input_hash] = M
    end
    M
end
