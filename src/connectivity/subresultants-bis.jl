# =========================================================================
# RING COERCION TOOLS (Replaces array wrangling and homogenize)
# =========================================================================

# TODO: use MPolyBuildCtx?

"""
    to_univariate(P::MPolyRingElem, var_idx::Int)

Converts a bivariate polynomial in `K[x, y]` into a univariate polynomial in `(K[x])[y]`.
"""
function to_univariate(P::MPolyRingElem, var_idx::Int)
    R = parent(P)
    K = base_ring(R)
    coeff_idx = 3 - var_idx # The variable to be pushed to the coefficients

    Rx, x = polynomial_ring(K, symbols(R)[coeff_idx])
    Rxy, y = polynomial_ring(Rx, symbols(R)[var_idx])

    res = zero(Rxy)
    for (c, exp) in zip(coefficients(P), exponent_vectors(P))
        # Build the coefficient in Rx, then attach to y in Rxy
        res += (c * x^exp[coeff_idx]) * y^exp[var_idx]
    end
    return res
end

"""
    to_bivariate(P::PolyRingElem, R_orig::MPolyRing)

Maps a univariate polynomial in `(K[x])[y]` back to the original bivariate ring `K[x, y]`.
"""
function to_bivariate(P::PolyRingElem, R_orig::MPolyRing)
    ctx = MPolyBuildCtx(R_orig)
    for (deg_y, c_y) in enumerate(coefficients(P))
        for (deg_x, c_x) in enumerate(coefficients(c_y))
            push_term!(ctx, c_x, [deg_x,deg_y])
        end
    end
    return finish(ctx)
end


# =========================================================================
# MODULAR ARITHMETIC & INTERSECTION
# =========================================================================

"""
    num_biv_rat_mod(A::MPolyRingElem, P::Vector)

Computes the numerator of `A(x, a(x)/b(x)) mod q(x)` efficiently using Residue Rings.
`P = [q, a, b]`.
"""
function num_biv_rat_mod(A::MPolyRingElem, P::Vector)
    q, a, b = P[1], P[2], P[3]
    Rx = parent(q)

    # Do math natively inside the quotient ring Rx / <q>
    U, mod_q = residue_ring(Rx, q)
    a_u, b_u, x_u = mod_q(a), mod_q(b), mod_q(gens(Rx)[1])

    deg_y = maximum(e[2] for e in exponent_vectors(A); init=0)
    res_u = zero(U)

    # Avoids homogenization entirely. b^(deg_y - j) clears the denominator.
    for (c, exp) in zip(coefficients(A), exponent_vectors(A))
        i, j = exp[1], exp[2]
        res_u += c * (x_u^i) * (a_u^j) * (b_u^(deg_y - j))
    end

    return lift(res_u)
end

"""
    intersect_biv(P::Vector, A::MPolyRingElem)

P = [q, a, b] encodes a finite set (x,a(x)/b(x)) where q(x)=0
Computes the divisor of `q` that encodes the intersection with `A(x,y)=0`.
Uses a modular CRT loop.
"""
function intersect_biv(P::Vector, A::MPolyRingElem)
    iszero(A) && return P[1]

    RS = symbols(parent(A))
    B, x = polynomial_ring(QQ, RS[1])

    dA_prev, dA_final = Int[], Int[]
    pprod, p = ZZ(1), ZZ(1) << 60
    compt = 0

    while compt < 12
        p = next_prime(p)
        Fp = GF(p)

        # Map polynomials natively to Fp
        Pp = [map_coefficients(Fp, poly) for poly in P]
        Ap = map_coefficients(Fp, A)

        Apev = num_biv_rat_mod(Ap, Pp)
        dAp = gcd(Pp[1], Apev)

        # Lift coefficients back to integers for CRT
        dA_current = [lift(ZZ, c) for c in coefficients(dAp)]

        if compt > 0
            dA_current = [crt([d1, d2], [pprod, p]) for (d1, d2) in zip(dA_prev, dA_current)]
        end

        count(!iszero, dA_current) <= 1 && return one(B) # Trivial GCD

        pprod *= p
        try
            dA_final = [reconstruct(c, pprod) for c in dA_current]
            (compt > 0 && dA_final == dA_prev) && break
        catch
            # Rational reconstruction failed, try next prime
        end

        dA_prev = dA_final
        compt += 1
    end

    return B(dA_final)
end


# =========================================================================
# SUBRESULTANTS
# =========================================================================

# Univariate resultant
function subresultants(P::Union{PolyRingElem{T}, FqPolyRingElem}, Q::Union{PolyRingElem{T}, FqPolyRingElem}, is_Fq::Bool=false) where T <: RingElem
    degree(P) < degree(Q) && ((P, Q) = (Q, P))

    S = typeof(P)[Q]
    s = leading_coefficient(Q)^(degree(P) - degree(Q)) # Initial scaling factor
    # rem is faster when working with finite fields (relies on FLINT)
    A = Q
    B = is_Fq ? (-leading_coefficient(Q))^(degree(P) - degree(Q) + 1) * rem(P, Q) : pseudorem(P, -Q)

    while !iszero(B)
        d, e = degree(A), degree(B)
        pushfirst!(S, B) # current remainder
        delta = d - e # degree drop

        if delta > 1
            if length(S) > 1
                n = degree(S[2]) - degree(S[1]) - 1 # stabilize scaling
                if n == 0
                    C = S[1]
                else
                    x, y = leading_coefficient(S[1]), leading_coefficient(S[2])
                    a = 1 << (length(bits(ZZ(n))) - 1)
                    c, n = x, n - a
                    while a > 1
                        a >>= 1
                        c = divexact(c^2, y)
                        if n >= a
                            c = divexact(c * x, y)
                            n -= a
                        end
                    end
                    C = divexact(c * S[1], y)
                end
            else
                # First step: fallback normalization
                C = divexact(leading_coefficient(B)^(delta - 1) * B, s^(delta - 1))
            end
            pushfirst!(S, C)
        else
            C = B # no degree drop
        end

        # Termination: constant polynomial
        e == 0 && return S

        # --- Next pseudo-remainder ---
        B = is_Fq ? (-leading_coefficient(B))^(d - e + 1) * rem(A, B) : pseudorem(A, -B)
        B = divexact(B, s^delta * leading_coefficient(A))

        A = C
        s = leading_coefficient(A)
    end
    return S
end

# Bivariate resultant
function subresultants(P::MPolyRingElem, Q::MPolyRingElem, var_idx::Int; list=false)
    UP, UQ = to_univariate(P, var_idx), to_univariate(Q, var_idx)
    sr_uni = subresultants(UP, UQ)

    R_orig = parent(P)
    #x_orig = symbols(R_orig)[3 - var_idx]
    v_orig = symbols(R_orig)

    if list
        # TODO: avoid evaluate, use change_ringvar?
        # Extract coefficients as bivariate polynomials evaluated at y=0
        return [change_ringvar(collect(coefficients(sr)), v_orig) for sr in sr_uni]
    else
        return [to_bivariate(sr, R_orig) for sr in sr_uni]
    end
end


# =========================================================================
# CORE ALGEBRAIC ENGINE
# =========================================================================

# delta : poly to factor w.r.t the polynomials in LP
function fact_gcd(delta::T, LP::Vector{T}) where T <: PolyRingElem
    Ldelta = Dict{Int, T}()
    curr_phi = gcd(delta, LP[1])
    Ldelta[1] = divexact(delta, curr_phi)

    for i in 2:length(LP)
        degree(curr_phi) == 0 && break
        next_phi = gcd(curr_phi, LP[i])
        Ldelta[i] = divexact(curr_phi, next_phi)
        curr_phi = next_phi
    end
    return filter(kv -> degree(kv.second) > 0, Ldelta)
end

# wrapper for the above function
function fact_gcd(delta::T, LP::Vector{T}) where (T <:MPolyRingElem)
    @assert is_univariate(delta) && all(is_univariate.(LP)) "Not univariate polynomial"
    R = parent(delta)
    K, RS = coefficient_ring(R), symbols(R)
    A, = polynomial_ring(K, RS[1])

    d_uni = A(coefficients_of_univariate(delta))
    LP_uni = [A(coefficients_of_univariate(p)) for p in LP]
    out = fact_gcd(d_uni, LP_uni)

    return Dict([ (i, change_ringvar(f, [first(RS)]) ) for (i,f) in out ])
end

"""
    param_crit_split(f::MPolyRingElem; v=0, detect_app=true)

Computes the critical points of a curve `f(x, y) = 0` grouped by multiplicity.
"""
function param_crit_split(f::MPolyRingElem, g::MPolyRingElem; v=0, detect_app=true)
    v > 0 && println("Compute subresultant sequence")
    f_y = derivative(f, 2)
    @iftime v>0 sr = subresultants(f, f_y, 2, list=true)

    (isempty(sr) || total_degree(sr[1][1]) == 0) && return Dict{Int, Vector}()

    v > 1 && println("Factorization")
    sqr_factors = factor(sr[1][1])

    # Filter non-zero multiplicity and sort
    sqr = sort([(fac, mult) for (fac, mult) in sqr_factors if mult > 0], by = x -> x[2])
    isempty(sqr) && return Dict{Int, Vector}()

    sqrmult = unique([s[2] for s in sqr])
    group_sqr = Dict(m => [s[1] for s in sqr if s[2] == m] for m in sqrmult)

    v > 0 && println("Compute crit partition w.r.t to multiplicity")
    param_crit = Dict{Int, Vector}()

    # Helper function to prevent out-of-bounds access on subresultant lists
    get_sr = (p, idx) -> (length(sr) >= p && length(sr[p]) >= abs(idx)) ? sr[p][end+idx] : zero(parent(f))

    singmult = filter(p -> p*(p-1) <= sqrmult[end], 2:sqrmult[end])
    for p in singmult
        param_crit[p] = [MPolyRingElem[], -get_sr(p, -1), (p-1)*get_sr(p, 0)]
    end

    # Critical points : multiplicity 1 in res
    (1 in sqrmult) && (param_crit[1] = [group_sqr[1], -get_sr(2, -1), get_sr(2, 0)])
    isempty(singmult) && return filter(p -> length(p.second[1]) > 0, param_crit)

    # Nodes : multiplicity 2 in res
    v > 0 && println("Compute apparent singularities")
    if 2 in sqrmult
        if detect_app
            f_x = derivative(f, 1)
            A = derivative(f_y, 2) * derivative(g, 1) - derivative(f_x, 2) * derivative(g, 2)

            dA = [intersect_biv([q, -get_sr(2, -1), get_sr(2, 0)], A) for q in group_sqr[2]]

            param_crit[0] = [group_sqr[2] ./ dA, -get_sr(2, -1), get_sr(2, 0)]
            append!(param_crit[2][1], dA)
        else
            param_crit[0] = [group_sqr[2], -get_sr(2, -1), get_sr(2, 0)]
        end
    end

    # Other singularities
    filter!(m -> !(m in [1, 2]), sqrmult)
    lsr = [get_sr(p, 0) for p in 2:(singmult[end]+1)]

    v > 0 && println("Compute gcd decomposition")
    for m in sqrmult
        for q in group_sqr[m]
            for (i, dji) in fact_gcd(q, lsr)
                haskey(param_crit, i+1) ? push!(param_crit[i+1][1], dji) : error("Curve not generic")
            end
        end
    end

    for v in values(param_crit)
        filter!(p -> total_degree(p) > 0, v[1])
    end

    return filter(p -> length(p.second[1]) > 0, param_crit)
end