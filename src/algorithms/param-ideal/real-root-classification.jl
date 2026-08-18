@doc Markdown.doc"""
    __pick_points(points::Vector{Vector{Vector{QQFieldElem}}}) -> Vector{Vector{QQFieldElem}}

Returns a list of sample points, one from each isolating box in `points`.

**Note**: This is an internal function.
"""
function _pick_points(points::Vector{Vector{Vector{QQFieldElem}}})::Vector{Vector{QQFieldElem}}
    return map(point -> map(x -> (x[1] + x[2]) / 2, point), points)
end

@doc Markdown.doc"""
    __unique_by_sign(points::Vector{Vector{QQFieldElem}}, f::Vector{QQMPolyRingElem}) -> Vector{Vector{QQFieldElem}}

Returns a list of points from `points` such that the sign vectors of `f` at these points are unique.

**Note**: This is an internal function.
"""
function _unique_by_sign(points::Vector{Vector{QQFieldElem}}, f::Vector{QQMPolyRingElem})::Vector{Vector{QQFieldElem}}
    signs = Set{Vector{QQFieldElem}}()
    res = Vector{Vector{QQFieldElem}}()
    for point in points
        s = [sign(evaluate(q, point)) for q in f]
        if !(s in signs)
            push!(signs, s)
            push!(res, point)
        end
    end
    res
end

@doc Markdown.doc"""
    _identity_matrix(m::Int) -> Matrix{Int}

Returns the identity matrix of size `m x m`.

**Note**: This is an internal function.
"""
function _identity_matrix(m::Int)::Matrix{Int}
    I = zeros(Int, m, m)
    for i in 1:m
        I[i, i] = 1
    end
    return I
end

@doc Markdown.doc"""
    _random_matrix(m::Int, n::Int) -> Matrix{Int}

Returns a random integer matrix of size `m x n` with entries in the range [1, 99].

**Note**: This is an internal function.
"""
function _random_matrix(m::Int, n::Int)::Matrix{Int}
    map(a -> 1 + abs(a) % 99, rand(Int, m, n))
end

@doc Markdown.doc"""
    real_root_classification(I::ParametricIdeal{K}, g::Vector{MPoly{K}}; <keyword arguments>) -> SemialgebraicSet

Computes the real root classification of the polynomials `g` with respect to the parametric ideal `I`. The result is a semi-algebraic set that describes the regions in the parameter space where the number of real roots such that `g` is positive is constant. Each region is accompanied by the constant number of real roots in that region, and a witness point in that region.

**Note**: The regions described by the semi-algebraic set may differ from the actual regions due to the fact that the Hermite matrices may not specialize well at certain parameter values. However, their symmetric difference is guaranteed to be contained in a proper algebraic set.

# Arguments
- `I::ParametricIdeal{K}`: input parametric ideal.
- `g::Vector{MPoly{K}}`: input list of polynomials.
- `nr_thrds::Int=1`: the number of threads to use for parallel computations.
- `worker_pool::AbstractWorkerPool=default_worker_pool()`: the worker pool to use for parallel computations.
- `info_level::Int=0`: info level printout: off (`0`, default), summary (`1`), detailed (`2`).
- `show_progress::Bool=info_level >= 1`: whether to show progress bars during computations.
- `output_form::Symbol=:sign`: the output form of the real root classification, can be either polynomials accompanied by their signs (`:sign`), Hermite matrices accompanied by their signatures (`:signature`), or Hermite matrices and multiplication matrices (`:matrix`).
- `ignore_no_real_roots::Bool=false`: whether to ignore sign conditions with no real roots.
"""
function real_root_classification(
    I::ParametricIdeal{K} where {K<:FracFieldElem},
    g::Vector{MPoly{K}} where {K<:FracFieldElem};
    nr_thrds::Int=1,
    worker_pool::AbstractWorkerPool=default_worker_pool(),
    info_level::Int=0,
    show_progress::Bool=info_level >= 1,
    output_form::Symbol=:sign,
    ignore_no_real_roots::Bool=false,
)::SemialgebraicSet
    R = base_ring(I.base_frac_field)
    # We first handle the case where a certain Hermite matrix is singular.
    # This can be handle by ideal saturation/radicalization.
    α = fill(0, length(g))
    H_1 = hermite_matrix(g, α, I; nr_thrds=nr_thrds, worker_pool=worker_pool, show_progress=show_progress)
    if size(H_1, 1) == 0
        return Empty()
    end
    h_1 = det(H_1)
    if is_zero(h_1)
        @warn "The Hermite matrix associated to one has zero determinant."
        @warn "Trying again with the radical ideal."
        I.radicalize = true
        return real_root_classification(I, g; nr_thrds=nr_thrds, worker_pool=worker_pool, show_progress=show_progress)
    end
    M = Vector{MatSpaceElem{FracFieldElem{QQMPolyRingElem}}}()
    m = Vector{FracFieldElem{QQMPolyRingElem}}()
    for i in 1:length(g)
        M_α = multiplication_matrix(I, g[i]; nr_thrds=nr_thrds, worker_pool=worker_pool, show_progress=show_progress)
        if M_α == zero_matrix(QQ, size(M_α, 1), size(M_α, 2)) && i == 1
            @warn "The multiplication matrix associated to the first polynomial is zero."
            @warn "This indicates an additional rank deficiency."
            @warn "Trying again by removing the first polynomial."
            return real_root_classification(I, g[2:end]; nr_thrds=nr_thrds, worker_pool=worker_pool, show_progress=show_progress)
        end
        m_α = det(M_α)
        if is_zero(m_α)
            @warn "The multiplication matrix associated to one polynomial has zero determinant."
            @warn "Trying again with the saturation ideal."
            push!(I.sats, g)
            return real_root_classification(I, g; nr_thrds=nr_thrds, worker_pool=worker_pool, show_progress=show_progress)
        end
        push!(M, M_α)
        push!(m, m_α)
    end
    if output_form == :matrix
        return MatEnum(hm=H_1, mm=M)
    end
    # We now compute a polynomial w whose vanishing defines a proper algebraic set of parameters,
    # where at least one of the Hermite matrices is singular or does not specialize well.
    # While this does not guarantee that the Hermite matrices specialize well outside of V(w),
    # the problematic parameters are contained in a proper algebraic set,
    # and thus the real root classification is generically correct.
    w = Vector{QQMPolyRingElem}()
    if is_empty(I.gens_alt)
        push!(w, _change_ring(denominator(h_1), R))
        push!(w, _change_ring(numerator(h_1), R))
        for m_α in m
            push!(w, _change_ring(denominator(m_α), R))
            push!(w, _change_ring(numerator(m_α), R))
        end
    else
        if info_level >= 2
            println("Computing GCD of determinants of Hermite matrices associated to the original and alternative generators.")
            println("Old factors: $(factor(numerator(h_1)))")
        end
        I.gens, I.gens_alt = I.gens_alt, I.gens
        H_1_alt = hermite_matrix(g, α, I; nr_thrds=nr_thrds, worker_pool=worker_pool, show_progress=show_progress)
        I.gens, I.gens_alt = I.gens_alt, I.gens
        h_1_alt = det(H_1_alt)
        push!(w, _change_ring(denominator(h_1), R))
        push!(w, _change_ring(denominator(h_1_alt), R))
        push!(w, _change_ring(gcd(numerator(h_1), numerator(h_1_alt)), R))
        if info_level >= 2
            println("New factors: $(factor(numerator(h_1_alt)))")
        end
        if info_level >= 1
            println("GCD of determinants of Hermite matrices: $(factor(w[end]))")
        end
        for m_α in m
            push!(w, _change_ring(denominator(m_α), R))
        end
    end
    # We now compute at least one point in each connected component of the semi-algebraic set
    # defined by the non-vanishing of the polynomials in w.
    # This allows us to determine the sign conditions of g at each connected component,
    # and thus to determine the adapted family of Hermite matrices that we need to consider.
    p = _pick_points(points_per_components(QQMPolyRingElem[], QQMPolyRingElem[], w; nr_thrds=nr_thrds, worker_pool=worker_pool, info_level=info_level))
    if info_level >= 1
        println("Number of sample points: ", length(p))
    end
    prog = Progress.ProgressBar(total=length(p); desc="Sign conditions", enabled=show_progress)
    Progress.update!(prog, 0)
    Σ = Vector{Vector{Vector{Int}}}([])
    c = Vector{Vector{Int}}([])
    for pᵢ in p
        Σᵢ, cᵢ = sign_determination(I, g, pᵢ)
        push!(Σ, Σᵢ)
        push!(c, cᵢ)
        Progress.next!(prog)
    end
    Progress.finish!(prog)
    if info_level >= 1
        println("Sign conditions at sample points: ", Set(Σ))
    end
    σ = fill(1, length(g))
    if !(σ in vcat(Σ...))
        if info_level >= 1
            println("Sign $σ not found at sample points; returning empty set.")
        end
        return Empty()
    end
    if info_level >= 1
        printstyled("[ Good news: ", bold=true, color=:light_magenta)
        println("sign $σ found at sample points!")
    end
    A = Vector{Vector{Int}}([])
    for i in 1:length(Σ)
        for j in i:length(Σ)
            append!(A, adapted_family(sort_signs(unique([Σ[i]; Σ[j]]))))
        end
    end
    A = sort_indices(unique(A))
    if info_level >= 1
        println("Adapted family for real root classification: ", A)
    end
    H = Vector{MatSpaceElem{FracFieldElem{QQMPolyRingElem}}}()
    for α in A
        H_α = hermite_matrix(g, α, I; nr_thrds=nr_thrds, worker_pool=worker_pool, show_progress=show_progress)
        push!(H, H_α)
    end
    if output_form == :signature
        sigs = Vector{Vector{Int}}()
        counts = Vector{Int}()
        witnesses = Vector{Vector{QQFieldElem}}()
        for (pᵢ, Σᵢ, cᵢ) in zip(p, Σ, c)
            j = findfirst(==(σ), Σᵢ)
            cᵢⱼ = isnothing(j) ? 0 : cᵢ[j]
            if ignore_no_real_roots && cᵢⱼ == 0
                continue
            end
            sig = [tarski_query(g, α, I, pᵢ; nr_thrds=nr_thrds, worker_pool=worker_pool, show_progress=show_progress) for α in A]
            if (!(sig in sigs))
                push!(sigs, sig)
                push!(counts, cᵢⱼ)
                push!(witnesses, pᵢ)
            end
        end
        return SigEnum(hms=H, sigs=sigs, counts=counts, witnesses=witnesses)
    end
    # We now compute the real root classification defined with semi-algebraic formulas.
    # The polynomials appearing in the semi-algebraic descriptions are leading principal minors
    # of the previously-identified Hermite matrices.
    # To determine their possible signs, another round of sample points is needed.
    Q = parent(H_1)(_identity_matrix(size(H_1, 1)))
    h = Vector{QQMPolyRingElem}()
    success = false
    while !success
        empty!(h)
        success = true
        for i in eachindex(A)
            α = A[i]
            if info_level >= 1
                println("Computing leading principal minors of Hermite matrix for: ", α)
            end
            H_α = transpose(Q) * H[i] * Q
            h_α = [det(H_α[1:j, 1:j]) for j in 1:size(H_α, 1)]
            if any(is_zero, h_α)
                @warn "A leading principal minor of the Hermite matrix for α = $(α) is zero."
                @warn "Applying a random change of basis to all Hermite matrices and trying again."
                # We could also try to generate Q deterministically,
                # but a random change of basis is simpler and works well in practice.
                Q = parent(H_1)(_random_matrix(size(H_1, 1), size(H_1, 1)))
                success = false
                break
            end
            append!(h, _change_ring(map(p -> numerator(p) * denominator(p), h_α), R))
        end
    end
    h = filter(!is_constant, h)
    p = _unique_by_sign(_pick_points(points_per_components(QQMPolyRingElem[], QQMPolyRingElem[], [w; h]; nr_thrds=nr_thrds, worker_pool=worker_pool, info_level=info_level)), h)
    if info_level >= 1
        println("Number of refined sample points: ", length(p))
    end
    prog = Progress.ProgressBar(total=length(p); desc="Refined sign conditions", enabled=show_progress)
    Progress.update!(prog, 0)
    signs = Vector{Vector{Int}}()
    counts = Vector{Int}()
    witnesses = Vector{Vector{QQFieldElem}}()
    for pᵢ in p
        Σᵢ, cᵢ = sign_determination(I, g, pᵢ)
        j = findfirst(==(σ), Σᵢ)
        cᵢⱼ = isnothing(j) ? 0 : cᵢ[j]
        if ignore_no_real_roots && cᵢⱼ == 0
            Progress.next!(prog)
            continue
        end
        new_sign = map(p -> Int(numerator(sign(evaluate(p, pᵢ)))), h)
        if !(new_sign in signs)
            push!(signs, new_sign)
            push!(counts, cᵢⱼ)
            push!(witnesses, pᵢ)
        end
        Progress.next!(prog)
    end
    Progress.finish!(prog)
    if info_level >= 1
        println("Signs at refined sample points: ", signs)
    end
    w = unique([q for wᵢ in w for (q, _) in factor(wᵢ)])
    # Finally, the semi-algebraic description of the real root classification is given by
    # the semi-algebraic set defined by the signs of the leading principal minors of the Hermite matrices,
    # restricted to the semi-algebraic set defined by the non-vanishing of the polynomials in w.
    Intersect([
        BasicSemialgebraicSet(eqs=[], ineqs=w, pos=[], nonneg=[]),
        SignEnum(polys=h, signs=signs, counts=counts, witnesses=witnesses)
    ])
end
