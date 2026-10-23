import IterTools: subsets

function incidence_variety(A::MatSpaceElem{<:MPolyRingElem}, r::Int, iota::Vector{Int})::Vector{<:MPolyRingElem}
    m = size(A, 1)
    @assert size(A) == (m, m)
    @assert length(iota) == m - r
    new_symbols = vec([Symbol("u_$(i)_$(j)") for j in 1:m-r, i in 1:m])
    R = base_ring(A)
    R′, vars = polynomial_ring(base_ring(R), [symbols(R); new_symbols], internal_ordering=internal_ordering(R))
    A′ = _change_ring(A, R′)
    new_vars = vars[number_of_variables(R)+1:end]
    Y = matrix_space(R′, m, m - r)(permutedims(reshape(new_vars, m - r, m), (2, 1)))
    return [vec(transpose(A′ * Y).entries); vec(transpose(Y[iota, :] - identity_matrix(QQ, m - r)).entries)]
end

function incidence_variety(A::MatSpaceElem{<:MPolyRingElem}, r::Int, U::Matrix{<:Real}, S::Matrix{<:Real})::Vector{<:MPolyRingElem}
    m = size(A, 1)
    @assert size(A) == (m, m)
    @assert size(U) == (m, m)
    @assert size(S) == (m, m)
    new_symbols = vec([Symbol("u_$(i)_$(j)") for j in 1:m-r, i in 1:m])
    R = base_ring(A)
    R′, vars = polynomial_ring(base_ring(R), [symbols(R); new_symbols], internal_ordering=internal_ordering(R))
    A′ = _change_ring(A, R′)
    U′ = matrix_space(R′, m - r, m)(U[1:m-r, 1:m])
    S′ = matrix_space(R′, m - r, m - r)(S[1:m-r, 1:m-r])
    new_vars = vars[number_of_variables(R)+1:end]
    Y = matrix_space(R′, m, m - r)(permutedims(reshape(new_vars, m - r, m), (2, 1)))
    return [vec(transpose(A′ * Y).entries); vec(transpose(U′ * Y - S′).entries)]
end

function jacobian(f::Vector{T}, t::Int=0)::MatSpaceElem{T} where {T<:MPolyRingElem}
    @assert !isempty(f)
    R = parent(f[1])
    r = length(f)
    n = number_of_variables(R) - t
    x = gens(R)[t+1:end]
    return matrix_space(R, r, n)(stack([[derivative(fᵢ, xⱼ) for fᵢ in f] for xⱼ in x]))
end

function lagrange_system(f::Vector{<:MPolyRingElem}, t::Int=0; use_minors::Bool=false)::Vector{<:MPolyRingElem}
    @assert !isempty(f)
    R = parent(f[1])
    p = length(f)
    if use_minors
        return [f; minors(jacobian(f, t)[:, 2:end], p)]
    end
    n = number_of_variables(R) - t
    new_symbols = [Symbol("l_$(i)") for i in 1:p]
    R′, vars = polynomial_ring(base_ring(R), [symbols(R); new_symbols], internal_ordering=internal_ordering(R))
    new_vars = vars[number_of_variables(R)+1:end]
    f′ = _change_ring(f, R′)
    z = matrix_space(R′, 1, p)(reshape(new_vars, 1, p))
    J = _change_ring(jacobian(f, t), R′)
    alpha = matrix_space(R′, 1, n)(reshape([one(R′); zeros(R′, n - 1)], 1, n))
    return [f′; vec((z * J - alpha).entries)]
end

function image(p::T, M::Matrix{<:Real}, t::Int=0)::T where {T<:MPolyRingElem}
    n = size(M, 1)
    @assert size(M) == (n, n)
    x = gens(parent(p))[t+1:t+n]
    evaluate(p, x, M * x)
end

function image(f::Vector{T}, M::Matrix{<:Real}, t::Int=0)::Vector{T} where {T<:MPolyRingElem}
    map(p -> image(p, M, t), f)
end

function fiber(f::Vector{T}, τ::Real, t::Int=0)::Vector{T} where {T<:MPolyRingElem}
    x = gens(parent(f[1]))[t+1]
    map(fᵢ -> evaluate(fᵢ, [x], [τ]), f)
end

function _lift(f::Vector{T}, τ::Real, t::Int=0)::Vector{T} where {T<:MPolyRingElem}
    x = gens(parent(f[1]))[t+1]
    return [x - τ; f]
end

function polar_varieties(f::Vector{U}, d::Int, M::Matrix{<:Real}, T::Vector{<:Real}, t::Int=0)::Vector{Vector{U}} where {U<:MPolyRingElem}
    @assert !isempty(f)
    if d < 0
        return Vector{Vector{U}}()
    end
    w₁ = [lagrange_system(image(f, M, t), t)]
    if d == 0
        return w₁
    end
    @assert !isempty(T)
    # Here, we play the trick of considering x₁ as a parameter,
    # so that the Jacobian is computed with respect to the remaining variables.
    w₂ = polar_varieties(fiber(f, T[1], t), d - 1, M[2:end, 2:end], T[2:end], t + 1)
    [w₁; map(w -> _lift(w, T[1], t), w₂)]
end

function expected_dimension(n::Int, m::Int, r::Int)::Int
    return n - (m - r + 1) * (m - r) // 2
end

function expected_minimum_rank(n::Int, m::Int)::Int
    for r in 0:m-1
        if expected_dimension(n, m, r) >= 0
            return r
        end
    end
    return m
end

function _reorder_variables(f::Vector{T}, n::Int)::Vector{T} where {T<:MPolyRingElem}
    if isempty(f)
        return f
    end
    R = parent(f[1])
    new_symbols = [symbols(R)[n+1:end]; symbols(R)[1:n]]
    R′, _ = polynomial_ring(base_ring(R), new_symbols, internal_ordering=internal_ordering(R))
    _change_ring(f, R′)
end

function real_det(A::MatSpaceElem{W}, M::Matrix{<:Real}, T::Vector{<:Real}, U::Matrix{<:Real}, S::Matrix{<:Real})::Vector{Vector{Vector{W}}} where {W<:MPolyRingElem}
    m = size(A, 1)
    @assert size(A) == (m, m)
    n = number_of_variables(base_ring(A))
    @assert size(M) == (n, n)
    @assert length(T) == n
    @assert size(U) == (m, m)
    @assert size(S) == (m, m)
    q = [Vector{Vector{W}}() for _ in 0:m-1]
    for r in 0:m-1
        d = expected_dimension(n, m, r)
        # To handle non-generic cases, you may set
        # d = n
        # However, please ensure that AlgebraicSolving.jl is not used
        # due to some bugs there!
        if d < 0
            continue
        end
        q[r+1] = polar_varieties(incidence_variety(A, r, U, S), d, M, T, 0)
        # To consider incidence varieties without generic choices of U and S,
        # uncomment the following lines and comment out the previous line.
        # for iota in subsets(1:m, m - r)
        #     append!(q[r+1], polar_varieties(incidence_variety(A, r, iota), d, M, T, t))
        # end
        q[r+1] = map(w -> _reorder_variables(w, n), q[r+1])
    end
    return q
end

function determinantal_variety(A::MatSpaceElem{T}, r::Int)::Vector{T} where {T<:MPolyRingElem}
    m = size(A, 1)
    @assert size(A) == (m, m)
    @assert 0 <= r < m
    return vec([det(A[rows, cols]) for rows in subsets(1:m, r + 1), cols in subsets(1:m, r + 1)])
end

function _generic_matrix(m::Int, n::Int)::Matrix{Int}
    M = zeros(Int, m, n)
    valid = false
    while !valid
        M = _random_matrix(m, n)
        valid = true
        for i in 1:min(m, n)
            if rank(M[1:i, 1:i]) < i || rank(M[(m-i+1):m, (n-i+1):n]) < i
                valid = false
                break
            end
        end
    end
    return M
end

function _random_point(n::Int)::Vector{Int}
    vec(_random_matrix(1, n))
end

function solve_lmi(
    A::MatSpaceElem{MPoly{FracFieldElem{QQMPolyRingElem}}};
    nr_thrds::Int=Threads.nthreads(),
    worker_pool::AbstractWorkerPool=default_worker_pool(),
    info_level::Int=2,
    compared::Bool=false
)
    m = size(A, 1)
    @assert size(A) == (m, m)
    @assert m > 0
    t = number_of_variables(base_ring(base_ring(base_ring(A))))
    n = number_of_variables(base_ring(A))
    println("Solving LMI with t = $t, n = $n, m = $m")
    println(repeat("-", 40))
    res = Vector{SemialgebraicSet}()
    # TODO: do we want to try _identity_matrix first?
    M = _generic_matrix(n, n)
    T = _random_point(n)
    U = _generic_matrix(m, m)
    S = _generic_matrix(m, m)
    g = image(alternate_signs(collect(coefficients(charpoly(A)))), M, 0)[1:(end-1)]
    g₀ = map(p -> constant_coefficient(evaluate(p, gens(base_ring(A))[1:n], fill(0, n))), g)
    g₀ = map(p -> denominator(p) * numerator(p), g₀)
    filter!(p -> !(is_constant(p) && sign(constant_coefficient(p)) >= 0), g₀)
    if isempty(g₀)
        printstyled("[ Good news: ", bold=true, color=:light_magenta)
        println("The LMI is feasible for all values of the parameters.")
        return Universe()
    end
    if !any(is_constant, g₀)
        push!(res, BasicSemialgebraicSet(eqs=[], ineqs=[], pos=[], nonneg=g₀))
    end
    f = real_det(A, M, T, U, S)
    f′ = if compared
        @warn "Computing a second family of polar varieties for comparison."
        M′ = _generic_matrix(n, n)
        T′ = _random_point(n)
        U′ = _generic_matrix(m, m)
        S′ = _generic_matrix(m, m)
        real_det(A, M′, T′, U′, S′)
    else
        nothing
    end
    for r in 0:(m-1)
        if isempty(f[r+1])
            continue
        end
        println("Solving for rank $r, with $(length(f[r+1])) polar varieties")
        for i in eachindex(f[r+1])
            fᵢ = f[r+1][i]
            fᵢ′ = compared ? f′[r+1][i] : Vector{MPoly{FracFieldElem{QQMPolyRingElem}}}()
            gᵢ = _change_ring(g, parent(fᵢ[1]))
            if r > 0
                d = _expected_degree(ParametricIdeal(fᵢ; num_vars=n))
                if d == 0
                    println("The $i-th polar variety has degree $d; skipping processing.")
                    continue
                end
                # sat = gᵢ[m-r+1]
                sat = sum(image(determinantal_variety(A, r - 1), M, 0) .^ 2)
                d_sat = _expected_degree(ParametricIdeal(fᵢ; num_vars=n, sats=[sat]))
                if d_sat == 0
                    println("The $i-th polar variety has degree $d; after removing low rank locus, degree becomes $d_sat; skipping processing.")
                    continue
                end
                # TODO: check if there is a rank deficiency
                println("The $i-th polar variety has degree $d; after removing low rank locus, degree becomes $d_sat")
                cond = real_root_classification(ParametricIdeal(fᵢ; num_vars=n, sats=[sat], gens_alt=fᵢ′), gᵢ[(m-r+1):end]; nr_thrds=nr_thrds, worker_pool=worker_pool, info_level=info_level)
                if cond isa Universe
                    printstyled("[ Good news: ", bold=true, color=:light_magenta)
                    println("The LMI is feasible for all values of the parameters.")
                    return Universe()
                end
                if !(cond isa Empty)
                    push!(res, cond)
                end
                println("Finished processing the $i-th polar variety.")
                println()
            else
                d = _expected_degree(ParametricIdeal(fᵢ; num_vars=n))
                if d == 0
                    println("The $i-th polar variety has degree $d; skipping processing.")
                    continue
                end
                println("The $i-th polar variety has degree $d")
                cond = real_root_classification(ParametricIdeal(fᵢ; num_vars=n, gens_alt=fᵢ′), gᵢ[(m-r+1):end]; nr_thrds=nr_thrds, worker_pool=worker_pool, info_level=info_level)
                if cond isa Universe
                    printstyled("[ Good news: ", bold=true, color=:light_magenta)
                    println("The LMI is feasible for all values of the parameters.")
                    return Universe()
                end
                if !(cond isa Empty)
                    push!(res, cond)
                end
                println("Finished processing the $i-th polar variety.")
                println()
            end
        end
        println(repeat("-", 40))
    end
    println("The LMI is now solved.")
    if isempty(res)
        return Empty()
    end
    UnionSet(res)
end

function solve_lmi(A::Matrix{MPoly{FracFieldElem{QQMPolyRingElem}}})
    m = size(A, 1)
    @assert size(A) == (m, m)
    @assert m > 0
    R = parent(A[1, 1])
    S = matrix_space(R, m, m)
    A′ = S(A)
    solve_lmi(A′)
end

function solve_lmi(A::Union{Matrix{Any},Matrix{MPolyRingElem}})
    R = nothing
    for i in axes(A, 1)
        for j in axes(A, 2)
            if isa(A[i, j], MPoly{FracFieldElem{QQMPolyRingElem}})
                if R === nothing
                    R = parent(A[i, j])
                elseif R != parent(A[i, j])
                    error("Entries of the input matrix belong to different rings.")
                end
            end
        end
    end
    if R === nothing
        error("No entry of the input matrix is a polynomial; cannot determine the base ring.")
    end
    A′ = map(a -> R(a), A)
    solve_lmi(A′)
end
