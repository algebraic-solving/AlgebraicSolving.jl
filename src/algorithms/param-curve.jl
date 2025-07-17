@doc Markdown.doc"""
    rational_curve_parametrization(I::Ideal{T} where T <: MPolyRingElem, <keyword arguments>)

Given an ideal `I` with solution set X being of dimension 1 over the complex numbers, return
a rational curve parametrization of the one-dimensional irreducible components of X.

In the output, the variables `x`,`y` of the parametrization correspond to the last two
entries of the `vars` attribute, in that order.

**Note**: At the moment only QQ is supported as ground field. If the dimension of the ideal
is not one an ErrorException is thrown.

# Arguments
- `I::Ideal{T} where T <: QQMPolyRingElem`: input generators.
- `info_level::Int=0`: info level printout: off (`0`, default), summary (`1`), detailed (`2`).
- `use_lfs::Bool=false`: add new variables (_Z2, _Z1) + 2 generic linear forms
- `cfs_lfs::Vector{Vector{ZZRingElem}} = []`: coefficients for the above linear forms
- `nr_thrds::Int=1`: number of threads for msolve
- `check_gen::Bool = true`: perform some genericity position checks on the last two variables

# Examples
```jldoctest
julia> using AlgebraicSolving

julia> R, (x1,x2,x3) = polynomial_ring(QQ, ["x1","x2","x3"])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x1, x2, x3])

julia> I = Ideal([x1+2*x2+2*x3-1, x1^2+2*x2^2+2*x3^2-x1])
QQMPolyRingElem[x1 + 2*x2 + 2*x3 - 1, x1^2 - x1 + 2*x2^2 + 2*x3^2]

julia> rational_curve_parametrization(I)
AlgebraicSolving.RationalCurveParametrization([:x1, :x2, :x3], Vector{ZZRingElem}[], x^2 + 4//3*x*y - 1//3*x + y^2 - 1//3*y, 4//3*x + 2*y - 1//3, QQMPolyRingElem[4//3*x^2 - 4//3*x*y + 2//3*x + 4//3*y - 1//3])

julia> rational_curve_parametrization(I, cfs_lfs=map.(Ref(ZZ),[[-8,2,2,-1,-8], [8,-7,-5,8,-7]]))
AlgebraicSolving.RationalCurveParametrization([:x1, :x2, :x3, :_Z2, :_Z1], Vector{ZZRingElem}[[-8, 2, 2, -1, -8], [8, -7, -5, 8, -7]], 4963//30508*x^2 - 6134//7627*x*y - 647//7627*x + y^2 + 1640//7627*y + 88//7627, -6134//7627*x + 2*y + 1640//7627, QQMPolyRingElem[8662//22881*x^2 - 21442//22881*x*y - 2014//7627*x + 9458//22881*y + 1016//22881, -2769//30508*x^2 + 4047//15254*x*y - 875//7627*x + 3224//7627*y + 344//7627, -9017//91524*x^2 + 9301//45762*x*y - 1185//7627*x + 8480//22881*y + 920//22881])
```
"""
function rational_curve_parametrization(
        I::Ideal{P} where P<:QQMPolyRingElem;                       # input generators
        info_level::Int=0,                                          # info level for print outs
        use_lfs::Bool = false,                                      # add generic variables
        cfs_lfs::Vector{Vector{ZZRingElem}} = Vector{ZZRingElem}[], # coeffs of linear forms
        nr_thrds::Int=1,                                            # number of threads (msolve)
        check_gen::Bool = true                                      # perform genericity check
    )
    info_level>0 && println("Compute ideal data" * (check_gen ? " and genericity check" : ""))
    lucky_prime = _generate_lucky_primes(I.gens, one(ZZ)<<30, one(ZZ)<<31-1, 1) |> first
    Itest = Ideal(change_base_ring.(Ref(GF(lucky_prime)), I.gens))
    Itest.dim = I.dim
    dimension(Itest)
    if Itest.dim == -1
        T = polynomial_ring(QQ, [:x,:y])[1]
        I.dim = -1
        I.rat_param = RationalCurveParametrization(Symbol[], Vector{ZZRingElem}[], T(-1), T(-1), QQMPolyRingElem[])
        return I.rat_param
    end
    @assert(Itest.dim==1, "I must define a curve or an empty set")
    if nvars(parent(I)) == 1
        T, (x,y) = polynomial_ring(QQ, [:x,:y])
        I.dim = 1
        I.rat_param = RationalCurveParametrization(Symbol[:x, :y], Vector{ZZRingElem}[[0,1]], y, T(1), QQMPolyRingElem[])
        return I.rat_param
    end

    DEG = hilbert_degree(Itest)
    # If generic variables must be added
    F, cfs_lfs = use_lfs ? _add_genvars(I.gens, max(2, length(cfs_lfs)), cfs_lfs) :
             !isempty(cfs_lfs) ? _add_genvars(I.gens, length(cfs_lfs), cfs_lfs) :
             (I.gens, Vector{ZZRingElem}[])
    R = parent(first(F))
    N = nvars(R)
    if check_gen let
        val = [ZZ(), ZZ()]
        # Bound on bifurcation set degree (e.g. Jelonek & Kurdyka, 2005)
        bif_bound = ZZ(1) << ( N * floor(Int, log2(maximum(total_degree.(F)))) + 1 )
        while any(is_divisible_by.(val, Ref(lucky_prime)))
            val = rand(-bif_bound:bif_bound, 2)
        end
        Fnew = vcat(F, val[1]*gens(R)[N-1] + val[2])
        new_lucky_prime = _generate_lucky_primes(Fnew, one(ZZ)<<30, one(ZZ)<<31-1, 1) |> first
        local INEW = Ideal(change_base_ring.(Ref(GF(new_lucky_prime)), Fnew))
        @assert(dimension(INEW) == 0 && hilbert_degree(INEW) == DEG, "The curve is not in generic position")
    end end

    # Compute DEG+2 evaluations of x in the param (whose total deg is bounded by DEG)
    PARAM  = Vector{Vector{QQPolyRingElem}}(undef, DEG+2)
    _values = Vector{ZZRingElem}(undef, DEG+2)
    i = 1
    free_ind = collect(1:DEG+2)
    used_ind = zeros(Bool, DEG+2)
    lc = nothing
    while length(free_ind) > 0
        if i > 2*(DEG+2)
            error("Too many bad specializations: permute variables or use_lfs=true")
        end
        # Evaluation of the generator at values x s.t. 0 <= |x|-i <= length(free_ind)/2
        # plus one point at -(length(free_ind)+1)/2 if the length if odd.
        # This reduces a bit the bitsize of the evaluation
        curr_values = ZZ.([-(i-1+(length(free_ind)+1)÷2):-i;i:(i-1+length(free_ind)÷2)])
        LFeval = Ideal.(_evalvar(F, N-1, curr_values))
        # Compute parametrization of each evaluation
        Lr = Vector{RationalParametrization}(undef, length(free_ind))
        for j in 1:length(free_ind)
            info_level>0 && print("Evaluated parametrizations: $(j+DEG+2-length(free_ind))/$(DEG+2)", "\r")
            Lr[j] = rational_parametrization(LFeval[j], nr_thrds=nr_thrds)
        end
        info_level>0 && println()
        for j in 1:length(free_ind)
            # Specialization checks: same vars order, generic degree
            if  Lr[j].vars == [symbols(R)[1:N-2]; symbols(R)[N]] && degree(Lr[j].elim) == DEG
                if isnothing(lc)
                    lc = leading_coefficient(Lr[j].elim)
                    rr = [ p for p in vcat(Lr[j].elim, Lr[j].denom, Lr[j].param) ]
                else
                    # Adjust when the rat_param is multiplied by some constant factor
                    fact = lc / leading_coefficient(Lr[j].elim)
                    rr = [ p*fact for p in vcat(Lr[j].elim, Lr[j].denom, Lr[j].param) ]
                end
                PARAM[j] = rr
                _values[j] = curr_values[j]
                used_ind[j] = true
            else
                info_level>0 && println("bad specialization: ", i+j-1)
            end
        end
        i += length(free_ind)
        free_ind = [ free_ind[j] for j in eachindex(free_ind) if !used_ind[j] ]
        used_ind = zeros(Bool, length(free_ind))
    end

    # Interpolate each coefficient of each poly in the param
    T = polynomial_ring(QQ, [:x,:y])[1]
    A = polynomial_ring(QQ)[1]

    POLY_PARAM = Vector{QQMPolyRingElem}(undef,N)
    for count in 1:N
        info_level>0 && print("Interpolate parametrizations: $count/$N\r")
        COEFFS = Vector{QQPolyRingElem}(undef, DEG+1)
        for deg in 0:DEG
            _evals = [coeff(PARAM[i][count], deg) for i in eachindex(PARAM)]
            # Remove denominators for faster interpolation with FLINT
            den = foldl(lcm, denominator.(_evals))
            scaled_evals = [ZZ(_evals[i] * den) for i in eachindex(_evals)]
            COEFFS[deg+1] = interpolate(A, _values, scaled_evals) / (lc*den)
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
    I.dim = 1 # If we reached here, I has necessarily dimension 1
    I.rat_param = RationalCurveParametrization( symbols(R), cfs_lfs, POLY_PARAM[1],
                                                POLY_PARAM[2], POLY_PARAM[3:end]    )
    return I.rat_param
end


# Inject polynomials in F in a polynomial ring with ngenvars new variables
# return these polynomials
# + newvars linear forms provided by coefficients in cfs_lfs + random ones
function _add_genvars(
    F::Vector{P} where P<:MPolyRingElem,
    ngenvars::Int,
    cfs_lfs::Vector{Vector{ZZRingElem}} = Vector{ZZRingElem}[]
)
    if length(cfs_lfs) > ngenvars
        error("Too many linear forms provided ($(length(cfs_lfs))>$(ngenvars))")
    end
    R = parent(first(F))
    N = nvars(R)
    # Add new variables (reverse alphabetical order)
    newS = vcat(symbols(R), Symbol.(["_Z$i" for i in ngenvars:-1:1]))
    R_ext, all_vars = polynomial_ring(base_ring(R), newS)
    # Inject F in this new ring
    Fnew = Vector{MPolyRingElem}(undef, length(F))
    ctx = MPolyBuildCtx(R_ext)
    for i in eachindex(F)
        for (e, c) in zip(exponent_vectors(F[i]), coefficients(F[i]))
            push_term!(ctx, c, vcat(e, zeros(Int,ngenvars)))
        end
        Fnew[i] = finish(ctx)
    end

    # Complete possible incomplete provided linear forms
    for i in eachindex(cfs_lfs)
        if N+ngenvars < length(cfs_lfs[i])
            error("Too many coeffs ($(length(cfs_lfs[i]))>$(N+ngenvars)) for the $(i)th linear form")
        else
            append!(cfs_lfs[i], rand(ZZ.(setdiff(-100:100,0)), N+ngenvars - length(cfs_lfs[i])))
        end
    end
    # Add missing linear forms if needed
    append!(cfs_lfs, [rand(ZZ.(setdiff(-100:100,0)), N+ngenvars) for _ in 1:ngenvars-length(cfs_lfs)])
    # Construct and append linear forms
    append!(Fnew, [transpose(cfs_lf) * all_vars for cfs_lf in cfs_lfs])

    return Fnew, cfs_lfs
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

# Generate N primes > start that do not divide any numerator/denominator
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