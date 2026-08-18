# Notation 2.43 (Sign variations)
@doc Markdown.doc"""
    sign_variations(s::Vector{T}) -> Int

Computes the number of sign variations in the sequence `s`, ignoring zeros.
"""
function sign_variations(s::Vector{T}) where {T}
    s′ = filter(x -> !is_zero(x), s)
    if isempty(s′)
        return 0
    end
    count = 0
    for i in 1:(length(s′)-1)
        if s′[i] * s′[i+1] < 0
            count += 1
        end
    end
    return count
end

@doc Markdown.doc"""
    alternate_signs(s::Vector{T}) -> Vector{T}

Returns a new vector where the signs of the elements in `s` are alternated, starting with the first element unchanged.
"""
function alternate_signs(s::Vector{T})::Vector{T} where {T}
    s′ = copy(s)
    for i in (length(s)%2+1):2:length(s)
        s′[i] = -s′[i]
    end
    return s′
end

@doc Markdown.doc"""
    signature(M::QQMatrix) -> Int

Computes the signature of the real symmetric matrix `M`, i.e., the number of positive eigenvalues minus the number of negative eigenvalues.
"""
function signature(M::QQMatrix)::Int
    @assert is_symmetric(M)
    l = collect(coefficients(charpoly(M)))
    sign_variations(l) - sign_variations(alternate_signs(l))
end

@doc Markdown.doc"""
    matrix_of_signs(A::Vector{Vector{Int}}, Σ::Vector{Vector{Int}}) -> Array{Int,2}

Computes the matrix of signs of `A` on `Σ`.
"""
function matrix_of_signs(A::Vector{Vector{Int}}, Σ::Vector{Vector{Int}})::Array{Int,2}
    m = length(A)
    n = length(Σ)
    M = ones(Int, m, n)
    for i in 1:m
        for j in 1:n
            for (αᵢ, σⱼ) in zip(A[i], Σ[j])
                if αᵢ == 0
                    continue
                elseif σⱼ == 0
                    M[i, j] = 0
                elseif αᵢ == 1
                    M[i, j] *= σⱼ
                elseif αᵢ == -1
                    M[i, j] /= σⱼ
                end
            end
        end
    end
    M
end

@doc Markdown.doc"""
    sort_indices(A::Vector{Vector{Int}}) -> Vector{Vector{Int}}

Sorts the list of multi-indices `A`.
"""
function sort_indices(A::Vector{Vector{Int}})::Vector{Vector{Int}}
    sort(A, by=α -> sum(α[i] * 3^(i - 1) for i in eachindex(α)))
end

@doc Markdown.doc"""
    sort_signs(Σ::Vector{Vector{Int}}) -> Vector{Vector{Int}}

Sorts the list of sign vectors `Σ`.
"""
function sort_signs(Σ::Vector{Vector{Int}})::Vector{Vector{Int}}
    sort(Σ, by=σ -> map(σᵢ -> σᵢ == -1 ? 2 : σᵢ, σ))
end

@doc Markdown.doc"""
    extend_signs(Σ::Vector{Vector{Int}}, T::Vector{Int}) -> Vector{Vector{Int}}

Extends the list of sign vectors `Σ` by appending each sign in `T` to each vector in `Σ`.
"""
function extend_signs(Σ::Vector{Vector{Int}}, T::Vector{Int})::Vector{Vector{Int}}
    Σ′ = Vector{Vector{Int}}()
    for σ in Σ
        for τ in T
            push!(Σ′, [σ; τ])
        end
    end
    Σ′
end

@doc Markdown.doc"""
    extend_indices(A::Vector{Vector{Int}}, B::Vector{Int}) -> Vector{Vector{Int}}

Extends the list of multi-indices `A` by appending each element in `B` to each vector in `A`.
"""
function extend_indices(A::Vector{Vector{Int}}, B::Vector{Int})::Vector{Vector{Int}}
    A′ = Vector{Vector{Int}}()
    for β in B
        for α in A
            push!(A′, [α; β])
        end
    end
    A′
end

@doc Markdown.doc"""
    modified_tensor_product(M::Matrix{Int}, M′::Matrix{Int}) -> Matrix{Int}

Computes the modified tensor product of the matrices `M` and `M′`.
"""
function modified_tensor_product(M::Matrix{Int}, M′::Matrix{Int})
    n, m = size(M)
    n′, m′ = size(M′)
    M′′ = zeros(Int, n * n′, m * m′)
    for i in 1:n
        for j in 1:m
            for i′ in 1:n′
                for j′ in 1:m′
                    M′′[n*(i′-1)+i, m′*(j-1)+j′] = M[i, j] * M′[i′, j′]
                end
            end
        end
    end
    M′′
end

@doc Markdown.doc"""
    adapted_family(Σ::Vector{Vector{Int}}) -> Vector{Vector{Int}}

Computes a family of multi-indices adapted to sign determination on the list of sign vectors `Σ`.
"""
function adapted_family(Σ::Vector{Vector{Int}})::Vector{Vector{Int}}
    if isempty(Σ)
        return []
    end
    s = length(Σ[1])
    if s == 0
        return [[]]
    end
    Ξ⁽⁰⁾ = Vector{Vector{Int}}([])
    Ξ⁽¹⁾ = Vector{Vector{Int}}([])
    Ξ⁽²⁾ = Vector{Vector{Int}}([])
    for σ in Σ
        τ = σ[1:(s-1)]
        if τ in Ξ⁽⁰⁾ && τ in Ξ⁽¹⁾
            Ξ⁽²⁾ = push!(Ξ⁽²⁾, τ)
        elseif τ in Ξ⁽⁰⁾
            Ξ⁽¹⁾ = push!(Ξ⁽¹⁾, τ)
        else
            Ξ⁽⁰⁾ = push!(Ξ⁽⁰⁾, τ)
        end
    end
    [
        extend_indices(adapted_family(Ξ⁽⁰⁾), [0]);
        extend_indices(adapted_family(Ξ⁽¹⁾), [1]);
        extend_indices(adapted_family(Ξ⁽²⁾), [2])
    ]
end

@doc Markdown.doc"""
    tarski_query(g::Vector{MPoly{K}}, α::Vector{Int}, I::ParametricIdeal{K}, vals::Vector{QQFieldElem}=QQFieldElem[]; nr_thrds::Int=1, worker_pool::AbstractWorkerPool=default_worker_pool(), show_progress::Bool=false) -> Int

Computes the Tarski query of `g^α` with respect to the parametric ideal `I` at the parameter values `vals`.
"""
function tarski_query(
    g::Vector{MPoly{K}},
    α::Vector{Int},
    I::ParametricIdeal{K},
    vals::Vector{QQFieldElem}=QQFieldElem[];
    nr_thrds::Int=1,
    worker_pool::AbstractWorkerPool=default_worker_pool(),
    show_progress::Bool=false
)::Int where {K<:FracFieldElem}
    @assert length(vals) == I.num_params
    return signature(hermite_matrix(g, α, I, vals; nr_thrds=nr_thrds, worker_pool=worker_pool, show_progress=show_progress))
end

@doc Markdown.doc"""
    tarski_query(g::Vector{MPoly{K}}, A::Vector{Vector{Int}}, I::ParametricIdeal{K}, vals::Vector{QQFieldElem}=QQFieldElem[]; nr_thrds::Int=1, worker_pool::AbstractWorkerPool=default_worker_pool(), show_progress::Bool=false) -> Vector{Vector{Int}}

Computes the Tarski query of `g^α` with respect to the parametric ideal `I` at the parameter values `vals` for each multi-index `α` in `A`.
"""
function tarski_query(
    g::Vector{MPoly{K}},
    A::Vector{Vector{Int}},
    I::ParametricIdeal{K},
    vals::Vector{QQFieldElem}=QQFieldElem[];
    nr_thrds::Int=1,
    worker_pool::AbstractWorkerPool=default_worker_pool(),
    show_progress::Bool=false
)::Vector{Vector{Int}} where {K<:FracFieldElem}
    [[tarski_query(g, α, I, vals; nr_thrds=nr_thrds, worker_pool=worker_pool, show_progress=show_progress)] for α in A]
end

@doc Markdown.doc"""
    sign_determination(I::ParametricIdeal{K}, g::Vector{MPoly{K}}, vals::Vector{QQFieldElem}=QQFieldElem[]; nr_thrds::Int=1, worker_pool::AbstractWorkerPool=default_worker_pool(), show_progress::Bool=false) -> Tuple{Vector{Vector{Int}}, Vector{Int}}

Computes the realizable sign conditions of the polynomials `g` with respect to the parametric ideal `I` at the parameter values `vals`. The result is a tuple containing a list of sign vectors and a list of number of real roots corresponding to each sign vector.
"""
function sign_determination(
    I::ParametricIdeal{K},
    g::Vector{MPoly{K}},
    vals::Vector{QQFieldElem}=QQFieldElem[];
    nr_thrds::Int=1,
    worker_pool::AbstractWorkerPool=default_worker_pool(),
    show_progress::Bool=false
)::Tuple{Vector{Vector{Int}},Vector{Int}} where {K<:FracFieldElem}
    s = length(g)
    r = tarski_query(g, fill(0, s), I, vals; nr_thrds=nr_thrds, worker_pool=worker_pool, show_progress=show_progress)
    if r == 0
        return (Vector{Vector{Int}}([]), Vector{Int}())
    end
    Σ = [Vector{Int}()]
    c = [r]
    A = [Vector{Int}()]
    for i in 1:length(g)
        Σ′ = extend_signs(Σ, [0, 1, -1])
        A′ = extend_indices(A, [0, 1, 2])
        M = matrix(QQ, matrix_of_signs(A′, Σ′))
        T = matrix(QQ, tarski_query(g[1:i], A′, I, vals; nr_thrds=nr_thrds, worker_pool=worker_pool, show_progress=show_progress))
        c′ = inv(M) * T
        @assert all(denominator(c′[j, 1]) == 1 for j in 1:length(Σ′))
        Σ = [Σ′[j] for j in 1:length(Σ′) if !is_zero(c′[j, 1])]
        c = [Int(numerator(c′[j, 1])) for j in 1:length(Σ′) if !is_zero(c′[j, 1])]
        @assert sum(c) == r
        A = adapted_family(Σ)
    end
    (Σ, c)
end
