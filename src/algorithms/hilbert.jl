@doc Markdown.doc"""
    hilbert_series(I::Ideal{T}) where T <: MPolyRingElem

Compute the Hilbert series of a given polynomial ideal `I`.

Based on: Anna M. Bigatti, Computation of Hilbert-Poincaré series,
Journal of Pure and Applied Algebra, 1997.

**Notes**:
* This requires a Gröbner basis of `I`, which is computed internally if not already known.
* Significantly faster when internal_ordering is :degrevlex.

# Examples
```jldoctest
julia> using AlgebraicSolving

julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"]);

julia> I = Ideal([x*y,x*z,y*z]);

julia> hilbert_series(I)
(3*t - 1)//(t - 1)
```
"""
function hilbert_series(I::Ideal{T}) where T <: MPolyRingElem

    gb = get!(I.gb, 0) do
        groebner_basis(I, complete_reduction = true)
    end
    lead_exps = Vector{Vector{Int}}(undef, length(gb))
    Threads.@threads for i in eachindex(gb)
        lead_exps[i] = _lead_exp_ord(gb[i], :degrevlex)
    end
    return _hilbert_series_mono(lead_exps)
end

@doc Markdown.doc"""
    hilbert_degree(I::Ideal{T}) where T <: MPolyRingElem

Compute the degree of a given polynomial ideal `I` by first computing its Hilbert series.

**Note**: This requires a Gröbner basis of `I`, which is computed internally if not already known.

# Examples
```jldoctest
julia> using AlgebraicSolving

julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"]);

julia> I = Ideal([x*y,x*z,y*z]);

julia> hilbert_degree(I)
3
```
"""
function hilbert_degree(I::Ideal{T}) where T <: MPolyRingElem

    !isnothing(I.deg) && return I.deg
    I.deg = numerator(hilbert_series(I))(1) |> abs
    return I.deg
end

@doc Markdown.doc"""
    hilbert_dimension(I::Ideal{T}) where T <: MPolyRingElem

Compute the Krull dimension of a given polynomial ideal `I` by first computing its Hilbert series.

**Note**: This requires a Gröbner basis of `I`, which is computed internally if not already known.

# Examples
```jldoctest
julia> using AlgebraicSolving

julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"]);

julia> I = Ideal([x*y,x*z,y*z]);

julia> hilbert_dimension(I)
1
```
"""
function hilbert_dimension(I::Ideal{T}) where T <: MPolyRingElem

    H = hilbert_series(I)
    I.dim = iszero(H) ? -1 : degree(denominator(H))
    return I.dim
end

@doc Markdown.doc"""
    hilbert_polynomial(I::Ideal{T}) where T <: MPolyRingElem

Compute the Hilbert polynomial and the index of regularity of a given polynomial ideal `I`
by first computing its Hilbert series. The index of regularity is the smallest integer such that
the Hilbert function and polynomial match.

Note that the Hilbert polynomial of I has leading term (e/d!)*t^d, where e and d are respectively
the degree and Krull dimension of I.

**Note**: This requires a Gröbner basis of `I`, which is computed internally if not already known.

# Examples
```jldoctest
julia> using AlgebraicSolving

julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"]);

julia> I = Ideal([x*y,x*z,y*z]);

julia> hilbert_degree(I)
(3*s + 3, 1)
```
"""
function hilbert_polynomial(I::Ideal{T}) where T <: MPolyRingElem

    A, s = polynomial_ring(QQ, :s)
    H = hilbert_series(I)
    dim = degree(denominator(H))
    num = iseven(dim) ? numerator(H) : -numerator(H)
    dim==0 && return num(s), 0

    t = gen(parent(num))
    La = Vector{ZZPolyRingElem}(undef, dim)
    while dim>0
        num, La[dim] = divrem(num, 1-t)
        dim -= 1
    end

    Hpolyfct = d->sum(La[i](0)*binomial(i+d, i) for i in 1:length(La))
    dim = degree(denominator(H))
    Hpoly = interpolate(A, QQ.(0:dim+1), [QQ(Hpolyfct(d)) for d in 0:dim+1])
    @assert(degree(Hpoly)==dim, "Degree of poly does not match the dimension")
    # Hilbert poly, index of regularity
    return Hpoly, degree(num)+1
end

# Computes hilbert series of a monomial ideal on input list of exponents
function _hilbert_series_mono(exps::Vector{Vector{Int}})

    h = _num_hilbert_series_mono(exps)
    t = gen(parent(h))
    return h//(1-t)^length(first(exps))
end

# Computes numerator hilbert series of a monomial ideal on input list of exponents
function _num_hilbert_series_mono(exps::Vector{Vector{Int}})

    A, t = polynomial_ring(ZZ, 't')
    r = length(exps)
    r == 0 && return one(A)
    N = length(first(exps))
    ## Base cases ##
    r == 1 && return (1-t^sum(first(exps)))
    supp = findall.(Ref(!iszero), exps)
    pow_supp = findall(s->length(s)==1, supp)
    # If exps is a product of simple powers
    if length(pow_supp) == r
        #println("Simple power")
        return prod(1-t^(exps[i][supp[i][1]]) for i in pow_supp)
    # Only one non-simple power P
    elseif length(pow_supp) == r-1
        #println("Mixed pow")
        inpow = setdiff(eachindex(exps), pow_supp) |> first
        # P has disjoint support with other powers
        if all(iszero(exps[inpow][supp[i][1]]) for i in pow_supp)
            return (1-t^sum(exps[inpow]))*prod(1-t^(exps[i][supp[i][1]]) for i in pow_supp)
        else
            return prod(1-t^(exps[i][supp[i][1]]) for i in pow_supp) - t^sum(exps[inpow]) *
            prod(1-t^(exps[i][supp[i][1]]-exps[inpow][supp[i][1]]) for i in pow_supp)
        end
    end

    # Variable index occuring the most in exps
    counts = sum(x->x .> 0, eachcol(reduce(hcat, exps)))
    ivarmax = argmax(counts)

    ## Splitting recursive cases ##
    # Monomials have disjoint supports
    if counts[ivarmax] == 1
        return prod(1-t^sum(mono) for mono in exps)
    # Heuristic where general splitting is useful
    elseif 8 <= r <= N
        # Finest partition of monomial supports
        LV, h = _monomial_support_partition(exps), one(A)
        rem_mon = collect(1:r)
        # If we are indeed splitting
        if length(LV) > 1
            for V in LV
                JV, iJV = Vector{Vector{Int}}(), Int[]
                for (k, i) in enumerate(rem_mon)
                    mono = exps[i]
                    if any(mono[j] != 0 for j in V)
                        push!(iJV, k)
                        push!(JV, mono)
                    end
                end
                h *= _num_hilbert_series_mono(JV)
                # Avoid re-check monomials
                deleteat!(rem_mon, iJV)
            end
            return h
        end
    end

    ## Pivot recursive case ##
    # Exponent of ivarmax in gcd of two random generators
    pivexp = max(1, minimum(mon[ivarmax] for mon in rand(exps, 2)))
    h = zero(A)
    #Compute and partition gens of (exps):pivot
    Lquo = [Vector{Int64}[] for _ in 1:pivexp+2]
    trivialquo = false
    for mono in exps
        if mono[ivarmax] <= pivexp
            monoquo = vcat(mono[1:ivarmax-1], 0, mono[ivarmax+1:end])
            if iszero(monoquo)
                trivialquo = true
                break
            end
            push!(Lquo[mono[ivarmax]+1],  monoquo)
        else
            push!(Lquo[pivexp+2],
                vcat(mono[1:ivarmax-1], mono[ivarmax]-pivexp, mono[ivarmax+1:end]))
        end
    end
    if !trivialquo
        # Interreduce generators based on partition
        @inbounds for i in pivexp+1:-1:1
            non_min = [ k for (k,mono) in enumerate(Lquo[i]) if
                    any(_all_lesseq(mini, mono) for j in i+1:pivexp+1 for mini in Lquo[j])]
            deleteat!(Lquo[i], non_min)
        end
        # Merge all partitions
        h += _num_hilbert_series_mono(vcat(Lquo...))*t^pivexp
    end
    # Interreduce (exps) + pivot
    filter!(e->(pivexp > e[ivarmax]), exps)
    push!(exps,[zeros(Int64,ivarmax-1); pivexp; zeros(Int64,N-ivarmax)])
    h += _num_hilbert_series_mono(exps)

    return h
end

function _all_lesseq(a::Vector{Int}, b::Vector{Int})::Bool
    @inbounds for i in eachindex(a)
        if a[i] > b[i]
            return false
        end
    end
    return true
end

# Build adjacency graph: connect variables that appear together in a monomial
function _monomial_support_partition(L::Vector{Vector{Int}})

    n = length(first(L))
    adj = [Set{Int}() for _ in 1:n]
    active = falses(n)

    for mono in L
        support = findall(!=(0), mono)
        foreach(i -> active[i] = true, support)
        for i in support, j in support
            i != j && push!(adj[i], j)
        end
    end

    visited = falses(n)
    components = Vector{Vector{Int}}()

    function dfs(u, comp)
        visited[u] = true
        push!(comp, u)
        foreach(v -> !visited[v] && dfs(v, comp), adj[u])
    end

    for v in 1:n
        if active[v] && !visited[v]
            comp = Int[]
            dfs(v, comp)
            push!(components, comp)
        end
    end

    return components
end
