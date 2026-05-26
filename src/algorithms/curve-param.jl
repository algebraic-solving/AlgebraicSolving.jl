@doc Markdown.doc"""
    curve_rational_parametrization(I::Ideal{T} where T <: MPolyRingElem, <keyword arguments>)

Given a **radical** ideal `I` with solution set X being of dimension 1 over the complex numbers,
return a rational curve parametrization of the one-dimensional irreducible components of X.

**Important**: In the output, the variables x and y correspond respectively to the last and second-to-last entries of the vars attribute.

**Note**: At the moment only QQ is supported as ground field. If the dimension of the ideal
is not one an ErrorException is thrown.

# Arguments
- `I::Ideal{T} where T <: QQMPolyRingElem`: input generators.
- `info_level::Int=0`: info level printout: off (`0`, default), summary (`1`), detailed (`2`).
- `cfs_lfs::Vector{Vector{ZZRingElem}} = []`: coefficients for the above linear forms
- `nr_thrds::Int=1`: number of threads for msolve

# Examples
```jldoctest
julia> using AlgebraicSolving

julia> R, (x1,x2,x3) = polynomial_ring(QQ, ["x1","x2","x3"])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x1, x2, x3])

julia> I = Ideal([x1+2*x2+2*x3-1, x1^2+2*x2^2+2*x3^2-x1])
QQMPolyRingElem[x1 + 2*x2 + 2*x3 - 1, x1^2 - x1 + 2*x2^2 + 2*x3^2]

julia> curve_rational_parametrization(I)
AlgebraicSolving.CurveRationalParametrization([:x1, :x2, :x3, :_Z2, :_Z1], Vector{ZZRingElem}[[0, 0, 1, 0, -1], [0, 1, 0, -1, 0]], x^2 + 4//3*x*y - 1//3*x + y^2 - 1//3*y, 4//3*x + 2*y - 1//3, QQMPolyRingElem[4//3*x^2 - 4//3*x*y + 2//3*x + 4//3*y - 1//3, -2*x^2 - 4//3*x*y + 2//3*x + 1//3*y, 4//3*x^2 + 2*x*y - 1//3*x])

julia> curve_rational_parametrization(I, cfs_lfs=map.(Ref(ZZ),[[-3,2,2], [1,4,-3]]))
AlgebraicSolving.CurveRationalParametrization([:x1, :x2, :x3, :_Z2, :_Z1], Vector{ZZRingElem}[[-3, 2, 2, 0, -1], [1, 4, -3, -1, 0]], 127//128*x^2 + 3//8*x*y + 161//64*x + y^2 - 7//8*y - 49//128, 3//8*x + 2*y - 7//8, QQMPolyRingElem[-3//32*x^2 - 1//2*x*y + 5//16*x + 1//2*y - 7//32, -1//4*x^2 + 1//8*x*y - 3//4*x + 3//8*y, 19//64*x^2 + 1//8*x*y + 25//32*x + 3//8*y - 21//64])

julia> curve_rational_parametrization(I, cfs_lfs=map.(Ref(ZZ),[[-3,2,2,-1,-2], [1,4,-3,2,-1]]))
AlgebraicSolving.CurveRationalParametrization([:x1, :x2, :x3, :_Z2, :_Z1], Vector{ZZRingElem}[[-3, 2, 2, -1, -2], [1, 4, -3, 2, -1]], 244//181*x^2 - 148//543*x*y + 532//543*x + y^2 + 182//181*y - 49//543, -148//543*x + 2*y + 182//181, QQMPolyRingElem[440//543*x^2 - 580//543*x*y - 44//543*x + 136//181*y + 112//543, 80//181*x^2 + 320//543*x*y + 122//181*x + 81//181*y + 49//543, -460//543*x^2 - 10//181*x*y - 418//543*x + 32//181*y + 56//181])
```
"""
function curve_rational_parametrization(
        I::Ideal{P} where P<:QQMPolyRingElem;                       # input generators
        info_level::Int=0,                                          # info level for print outs
        cfs_lfs::Vector{Vector{ZZRingElem}} = Vector{ZZRingElem}[], # coeffs of linear forms
        nr_thrds::Int=1,                                            # number of threads (msolve)
    )
    @assert(nvars(parent(I))>=2, "I must be defined in a ring with at least 2 variables")

    info_level>0 && println("Compute generic linear forms...")
    Inew, cfs_lfs = _add_genvars(I, 2, cfs_lfs)

    if Inew.dim == -1
        T = polynomial_ring(QQ, [:x,:y])[1]
        I.dim = -1
        I.rat_param = CurveRationalParametrization(Symbol[], Vector{ZZRingElem}[], T(-1), T(-1), QQMPolyRingElem[])
        return I.rat_param
    end
    @assert Inew.dim == 1 "I must define a curve or an empty set"

    R = parent(Inew)
    N = nvars(R)
    DEG, F = Inew.deg, Inew.gens

    # Compute DEG+2 evaluations of x in the param (whose total deg is bounded by DEG)
    PARAM  = Vector{Vector{QQPolyRingElem}}(undef, DEG+2)
    _values = Vector{ZZRingElem}(undef, DEG+2)

    i = 1
    free_ind = collect(1:DEG+2)
    used_ind = falses(DEG+2)
    lc = nothing

    while length(free_ind) > 0
        if i > 2*(DEG+2)
            error("Too many bad specializations. Check radicality and generic linear forms.")
        end
        # Determine values to evaluate at to keep bitsize low
        curr_values = ZZ.([-(i - 1 + (length(free_ind) + 1) ÷ 2):-i;i:(i - 1 + length(free_ind) ÷ 2)])
        LFeval = Ideal.(_evalvar(F, N, curr_values))
        # Compute parametrization of each evaluation
        Lr = Vector{RationalParametrization}(undef, length(free_ind))

        for j in eachindex(free_ind)
            info_level>0 && print("Evaluated parametrizations: $(j)/$(length(free_ind))", "\r")
            Lr[j] = rational_parametrization(LFeval[j], nr_thrds=nr_thrds)

            # Specialization checks: same vars order, generic degree
            if  Lr[j].vars == symbols(R)[1:N-1] && degree(Lr[j].elim) == DEG
                if isnothing(lc)
                    lc = leading_coefficient(Lr[j].elim)
                    rr = vcat([Lr[j].elim, Lr[j].denom], Lr[j].param)
                else
                    # Adjust when the rat_param is multiplied by some constant factor
                    fact = lc / leading_coefficient(Lr[j].elim)
                    rr = vcat([Lr[j].elim * fact, Lr[j].denom * fact], Lr[j].param .* fact)
                end
                PARAM[free_ind[j]] = rr
                _values[free_ind[j]] = curr_values[j]
                used_ind[j] = true
            end
        end

        # Update range, free indices and used indices
        i += length(free_ind)
        free_ind = free_ind[.!used_ind]
        used_ind = falses(length(free_ind))

        info_level * length(free_ind) != 0 &&
        println("bad specialization(s): ", curr_values[free_ind])
    end

    # Interpolate each coefficient of each poly in the param
    T, = polynomial_ring(QQ, [:x,:y])
    A, = polynomial_ring(QQ)

    POLY_PARAM = Vector{QQMPolyRingElem}(undef,N)
    for count in 1:N
        info_level>0 && print("Interpolate parametrizations: $count/$N\r")
        COEFFS = Vector{QQPolyRingElem}(undef, DEG+1)
        for deg in 0:DEG
            _evals = [coeff(PARAM[i][count], deg) for i in eachindex(PARAM)]
            # Remove denominators for faster interpolation with FLINT
            # TODO: remove dens mult when interface's ready in Nemo
            den = foldl(lcm, denominator.(_evals))
            scaled_evals = [ZZ(_evals[i] * den) for i in eachindex(_evals)]
            COEFFS[deg+1] = interpolate(A, _values, scaled_evals) / (lc * den)
        end

        ctx = MPolyBuildCtx(T)
        for (i, c) in enumerate(COEFFS)
            for (j, coeff) in enumerate(coefficients(c))
                !iszero(coeff) && push_term!(ctx, coeff, [j-1, i-1])
            end
        end
        POLY_PARAM[count] = finish(ctx)
    end
    info_level>0 && println()

    # Output: [vars, linear forms, elim, denom, [nums_param]]
    I.deg, I.dim = Inew.deg, Inew.dim
    I.rat_param = CurveRationalParametrization( symbols(R), cfs_lfs, POLY_PARAM[1],
                                                POLY_PARAM[2], POLY_PARAM[3:end]    )
    return I.rat_param
end


# Return F in a polynomial ring with ngenvars new variables
# + newvars linear forms provided by coefficients in cfs_lfs or generic ones internally computed
function _add_genvars(
    I::Ideal{P} where P<:MPolyRingElem,
    ngenvars::Int,
    cfs_lfs::Vector{Vector{T}} where T<:RingElem
)
    F = I.gens
    R = parent(I)
    K, n = base_ring(R), nvars(R)

    # Add new variables (reverse index order)
    newS = vcat(symbols(R), Symbol.(["_Z$i" for i in ngenvars:-1:1]))
    R_ext, all_vars = polynomial_ring(K, newS)

    # Inject F in this new ring efficiently using evaluation
    Fnew = Vector{MPolyRingElem}(undef, length(F))
    ctx = MPolyBuildCtx(R_ext)
    for i in eachindex(F)
        for (e, c) in zip(exponent_vectors(F[i]), coefficients(F[i]))
            push_term!(ctx, c, vcat(e, zeros(Int, ngenvars)))
        end
        Fnew[i] = finish(ctx)
    end

    # Find generic linear forms
    if !isempty(cfs_lfs)
        @assert length(cfs_lfs) == ngenvars "Expected $ngenvars linear forms, got $(length(cfs_lfs))"
        @assert all(length(c) in [n, n + ngenvars] for c in cfs_lfs) "Linear forms must have $n or $(n + ngenvars) coefficients"
    end
    if characteristic(K) == 0
        (DEG, DIM), cfs_lfs = _find_generic_linear_forms(I, ngenvars, cfs_lfs)
    else
        DEG, DIM = hilbert_degree(I), dimension(I)
    end

    # Extend coefficient vectors if only X-coefficients are provided
    if length(first(cfs_lfs)) == n
        cfs_lfs = [
            vcat(c, [-ZZ(j == i) for j in ngenvars:-1:1])
            for (i, c) in enumerate(cfs_lfs)
        ]
    end

    # Add equations Li(X) - Zi = 0
    append!(Fnew, [transpose(c) * all_vars for c in cfs_lfs])

    Inew = Ideal(Fnew)
    Inew.deg, Inew.dim = DEG, DIM

    return Inew, cfs_lfs
end

# Computes ngenvars sequential generic linear forms.
function _find_generic_linear_forms(I::Ideal{T}, ngenvars::Int, cfs_lfs) where T <: MPolyRingElem
    R, F = parent(I), I.gens
    n, vars = nvars(R), gens(R)
    cfs_lfs_out = Vector{Vector{ZZRingElem}}()

    # 1. Compute the degree of the system
    lucky_prime = first(_generate_lucky_primes(F, one(ZZ)<<30, (one(ZZ)<<31)-1, 1))
    if haskey(I.gb, 0)
        DEG, DIM = hilbert_degree(I), dimension(I)
    else
        Itest = Ideal(change_base_ring.(Ref(GF(lucky_prime)), F))
        DEG, DIM = hilbert_degree(Itest), dimension(Itest)
    end
    @assert DIM < 0 || ngenvars < DIM + 2 "Too many generic linear forms asked > dim + 1"

    !isempty(cfs_lfs) && return (DEG, DIM), cfs_lfs

    # 2. Compute a bound for generic specialization values
    # Bound on bifurcation set degree (e.g., Jelonek & Kurdyka, 2005)
    max_deg = maximum(f -> total_degree(f), F; init=1)
    bif_bound = ZZ(1) << (n * floor(Int, log2(max_deg)) + 1)

    # 3. Instantiate the stateful generator ONCE for all iterations
    candidate_stream = _candidate_stream(n)

    # Running system to test subsequent linear forms
    current_F = copy(F)
    for k in 1:ngenvars
        # 3. Compute a generic specialization value
        val = [ZZ(), ZZ()]
        while iszero(val[1]) || is_divisible_by(val[1], lucky_prime) || is_divisible_by(val[2], lucky_prime)
            val = rand(-bif_bound:bif_bound, 2)
        end

        # 4. Find the next generic linear form
        coeffs = _search_single_linear_form(current_F, val, DEG, DIM, k, candidate_stream)
        push!(cfs_lfs_out, coeffs)

        # 5. Specialize the current linear form
        L = sum(coeffs[i] * vars[i] for i in 1:n)
        push!(current_F, val[1] * L + val[2])

    end

    return (DEG, DIM), cfs_lfs_out
end


# Returns a generic linear form provided the current situation
# It gets the next candidate from the stream and applies genericity tests
function _search_single_linear_form(F, val, DEG, DIM, k, candidate_stream; max_iter=10000)
    R = parent(first(F))
    n, vars = nvars(R), gens(R)

   # Take candidates from the stream until we find a match
    for _ in 1:max_iter
        coeffs = take!(candidate_stream)

        L = sum(coeffs[i] * vars[i] for i in 1:n)
        Feval = vcat(F, val[1] * L + val[2])
        lucky_prime = first(_generate_lucky_primes(Feval, one(ZZ)<<30, (one(ZZ)<<31)-1, 1))
        Imod = Ideal(change_base_ring.(Ref(GF(lucky_prime)), k <= DIM ? Feval : F))

        if k <= DIM
            # --- Case A: Projection Form ---
            dimension(Imod) == DIM - k && hilbert_degree(Imod) == DEG && return coeffs
        else
            # --- Case B: Separating Form (k == DIM + 1) ---
            Iext, _ = _add_genvars(Imod, 1, [GF(lucky_prime).(coeffs)])
            Iext_elim = Ideal(eliminate(Iext, n))

            hilbert_degree(Iext_elim) == DEG && return coeffs
        end
    end

    error("Failed to find a generic linear form after $max_iter tests.")
end

# A stateful, lazy generator for candidate linear forms with n vars
function _candidate_stream(n::Int)
    Channel{Vector{ZZRingElem}}() do ch
        # 1. Quick coordinate projection check backward from the last variable
        for i in n:-1:1
            coeffs = zeros(ZZRingElem, n)
            coeffs[i] = 1
            put!(ch, coeffs)
        end

        # 2. Bitmask Layering Search
        sorted_masks = sort(collect(1:(1 << (n-1)) - 1), by=count_ones)
        queue = [ones(ZZRingElem, n)]
        tested = Set{Vector{ZZRingElem}}([queue[1]])

        while true
            # Yield the current layer's candidates one by one
            for coeffs in queue
                put!(ch, coeffs)
            end

            # Generate the next layer
            new_L = [
                l .+ [ZZ((mask >> (k-1)) & 1) for k in 1:n]
                for mask in sorted_masks
                for l in queue
            ]

            queue = Vector{Vector{ZZRingElem}}()
            for coeffs in new_L
                if !(coeffs in tested)
                    push!(tested, coeffs)
                    push!(queue, coeffs)
                end
            end
        end
    end
end

# for each a in La, evaluate each poly in F in x_i=a
function _evalvar(
    F::Vector{P} where P<:MPolyRingElem,
    i::Int,
    La::Vector{T} where T<:RingElem
    )
    R = parent(first(F))
    indnewvars = setdiff(1:nvars(R), i)
    C, = polynomial_ring(base_ring(R), symbols(R)[indnewvars])

    LFeval = Vector{typeof(zero(C))}[]
    ctx = MPolyBuildCtx(C)

    for a in La
        powa = Dict(0=>one(parent(a)), 1=>a) #no recompute powers
        push!(LFeval, typeof(zero(C))[])
        for f in F
            for (e,c) in zip(exponent_vectors(f), coefficients(f))
                aei = get!(powa, e[i]) do
                    a^e[i]
                end
                push_term!(ctx, c*aei, [e[j] for j in indnewvars ])
            end
            push!(LFeval[end], finish(ctx))
        end
    end
    return LFeval
end

# Generate N random primes between low and up
# that do not divide any numerator/denominator
# of any coefficient in polynomials from LP
function _generate_lucky_primes(
    LF::Vector{P} where P<:MPolyRingElem,
    low::ZZRingElem,
    up::ZZRingElem,
    N::Int64
    )
    # Avoid repetitive enumeration and redundant divisibility check
    CF = ZZRingElem[]
    for f in LF, c in coefficients(f), part in (numerator(c), denominator(c))
        if !isone(part)
            push!(CF, part)
        end
    end
    sort!(CF, rev=true)
    unique!(CF)

    # Test primes
    Lprim = ZZRingElem[]
    while length(Lprim) < N
        cur_prim = next_prime(rand(low:up))
        is_lucky = !(cur_prim in Lprim)
        i = firstindex(CF)
        # Exploit decreasing order of CF
        while is_lucky && i <= lastindex(CF) && CF[i] > cur_prim
            is_lucky = !is_divisible_by(CF[i], cur_prim)
            i += 1
        end
        is_lucky && push!(Lprim, cur_prim)
    end
    return Lprim
end
