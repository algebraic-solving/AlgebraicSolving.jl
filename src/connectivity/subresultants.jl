# Some tools

function homogenize(P::MPolyRingElem, i::Int, S::Vector{T} where T<:Union{Symbol, String})
    A = base_ring(parent(P))
    if iszero(P)
        return zero(polynomial_ring(A, S)[1])
    end
    # Homogenize w.r.t to the i-th variables
    LPc, LPe = poly_to_array(P)
    deg, N = maximum(getindex.(LPe, i)), length(LPe[1])
    for j in eachindex(LPe)
        push!(LPe[j], deg-LPe[j][i])
    end
    R, = polynomial_ring(A, S)
    return array_to_poly([LPc,LPe], R)
end

function homogenize(LP::Vector{P} where P<:MPolyRingElem, i::Int, S::Vector{T} where T<:Union{Symbol, String})
    return [ homogenize(lp, 2, S) for lp in LP ]
end

function poly_to_array(P::MPolyRingElem)
    # return P as two lists ci and vi of resp. coeffs and exponents
    return [collect(coefficients(P)), collect(exponent_vectors(P))]
end

function rem_ind(L,i)
    # remove element at index i
    return [ L[1:i-1]; L[i+1:end] ]
end

function add_ind(L,i,x)
    # add element x at index i
    return [ L[1:i-1]; [x]; L[i:end] ]
end

function parray_asvarcoeff(LP, idx)
    # takes the above representation of a poly P
    # and outputs a representation of P with univariate coeffs in the i-th variable
    A = parent(LP[1][1])
    NLP = [ Vector{typeof(LP[1][1])}[], [] ]
    for (c,e) in zip(LP[1], LP[2])
        vcat(zeros(A, e[idx]), [c])
        push!(NLP[1], vcat(zeros(A, e[idx]), [c]))
        push!(NLP[2], rem_ind(e, idx))
    end
    return NLP
end

function parray_asvar(LP, idx)
    # takes the above representation of a poly P
    # and outputs a representation of P seen as a univariate poly in the i-th variable
    deg = sort(LP[2], by= x->x[idx])[end][idx]
    NLP = [ [[],[]] for _ in 1:deg+1 ]
    for i in eachindex(LP[2])
        di = LP[2][i][idx]+1
        push!(NLP[di][1], LP[1][i])
        push!(NLP[di][2], rem_ind(LP[2][i], idx))
    end
    return NLP
end

function array_to_poly(LP, A)
    C = MPolyBuildCtx(A)
    R = base_ring(A)
    for i in eachindex(LP[1])
        push_term!(C, R(LP[1][i]), LP[2][i]);
    end
    return finish(C)
end

"""
    subresultants(P, Q)
Compute the subresultant polynomial remainder sequence (PRS) of two univariate polynomials.

Returns a vector `[S0, S1, ..., Sk]` where:
- S0 = Q (after possible swap)
- each Si is a subresultant
- sequence ends when remainder becomes zero or constant
"""

function subresultants(P::PolyRingElem{T}, Q::PolyRingElem{T}) where T <: RingElement
    return _subresultants(P, Q, false)
end

function subresultants(P::FqPolyRingElem, Q::FqPolyRingElem)
    return _subresultants(P, Q, true)
end

function _subresultants(P::Union{PolyRingElem{T}, FqPolyRingElem}, Q::Union{PolyRingElem{T}, FqPolyRingElem}, Fq::Bool) where T <: RingElement

    # Ensure deg(P) >= deg(Q)
    if degree(P) < degree(Q)
        P, Q = Q, P
    end

    S = Vector{typeof(Q)}()
    push!(S, Q)

    # Initial scaling factor (controls normalization of pseudo-remainders)
    s = leading_coefficient(Q)^(degree(P) - degree(Q))

    # rem is faster when working with finite fields (relies on FLINT)
    A = Q
    B = Fq ? (-leading_coefficient(Q))^(degree(P) - degree(Q) + 1) * rem(P, Q) : pseudorem(P, -Q)

    while true
        d = degree(A)
        e = degree(B)
        # Termination: exact division reached
        if iszero(B)
            return S
        end
        pushfirst!(S, B)  # prepend current remainder
        delta = d - e  # degree drop
        # --- Compute normalization factor C ---
        if delta > 1
            if length(S) > 1
                # Use previous subresultants to stabilize scaling
                n = degree(S[2]) - degree(S[1]) - 1
                if n == 0
                    C = S[1]
                else
                    x = leading_coefficient(S[1])
                    y = leading_coefficient(S[2])

                    # fast exponentiation (binary powering)
                    a = 1 << (length(bits(ZZ(n))) - 1)
                    c = x
                    n -= a
                    while a > 1
                        a >>= 1
                        c = c^2 / y
                        if n >= a
                            c = c * x / y
                            n -= a
                        end
                    end

                    C = c * S[1] / y
                end
            else
                # First step: fallback normalization
                C = leading_coefficient(B)^(delta - 1) * B / s^(delta - 1)
            end
            pushfirst!(S, C)
        else
            C = B
        end

        # Termination: constant polynomial
        if e == 0
            return S
        end

        # --- Next pseudo-remainder ---
        B = Fq ? (-leading_coefficient(B))^(d - e + 1)* rem(A, B) : pseudorem(A, -B)
        B /= (s^delta * leading_coefficient(A))

        # Shift sequence
        A = C
        s = leading_coefficient(A)
    end
end

#Bivariate subresultants
function subresultants(P::MPolyRingElem{T}, Q::MPolyRingElem{T}, idx; list=false) where T <: RingElement
    LPQ = map(poly_to_array, [P,Q])
    ULPQ = [ parray_asvar(lpq, idx) for lpq in LPQ ]

    R, x = polynomial_ring(base_ring(parent(P)), "x")
    S, y = polynomial_ring(R, "y")

    UP, UQ = [ S([ array_to_poly(l, R) for l in lpq ]) for lpq in ULPQ ]

    sr = subresultants(UP, UQ)

    # Get it back to initial polynomial ring
    Lsr = [ [collect(coefficients(csr)) for csr in coefficients(sr)] for sr in sr]
    newsr = []
    if list
        for (k,lsr) in enumerate(Lsr)
            mlsr = []
            for i in 1:k
                tmp = [[],[]]
                if i in eachindex(lsr)
                    for j in 1:length(lsr[i])
                        push!(tmp[1], lsr[i][j])
                        push!(tmp[2],add_ind([j-1],idx, 0))
                    end
                else
                    tmp = [[P|>parent|>base_ring|>zero], [[0,0]]]
                end
                push!(mlsr, tmp)
            end
            mlsr = [ array_to_poly(tmp, parent(P)) for tmp in mlsr ]
            push!(newsr, mlsr)
        end
    else
        for lsr in Lsr
            mlsr = [[],[]]
            for i in 1:length(lsr)
                for j in 1:length(lsr[i])
                    push!(mlsr[1], lsr[i][j])
                    push!(mlsr[2],add_ind([j-1],idx, i-1))
                end
            end
            push!(newsr, array_to_poly(mlsr, parent(P)))
        end
    end

    return newsr
end

function fact_gcd(delta::T, LP::Vector{T}) where (T <:PolyRingElem)
    # delta : poly to factor w.r.t the polynomials in LP
    Lphi = [gcd(delta, LP[1])]
    Ldelta = [delta/Lphi[1]]
    i = 2
    while degree(Lphi[end])>0
        #println.([Lphi,Ldelta,""])
        push!(Lphi, gcd(Lphi[i-1], LP[i]))
        push!(Ldelta, Lphi[i-1]/Lphi[i])
        i+=1
    end
    return filter(s->degree(s[2])>0, Dict(enumerate(Ldelta)))
end

function num_biv_rat_mod(A, P)
    # Computes the numerator of the polynomial A(x,a(x)//b(x)) mod q
    # where P = [q, a, b]
    B = base_ring(parent(A))
    RS = symbols(B)

    T, = polynomial_ring(B, RS[1])
    Puniv = [ T(coefficients_of_univariate(p)) for p in P ]

    U, fU = Nemo.residue_ring(T, Puniv[1])
    amod, bmod = fU.(Puniv[2:end])

    newv = vcat(RS, [:t])
    Tbiv, = polynomial_ring(T, newv[2:end])
    Ah = homogenize(A, 2, newv)
    Ah1 = array_to_poly(parray_asvarcoeff(poly_to_array(Ah), 1), Tbiv)
    Ahmod = change_coefficient_ring(U, Ah1)

    Aeval = evaluate(Ahmod, [amod, bmod])
    return change_ringvar(lift(Aeval), RS)
end


function intersect_biv(P::Vector{T} where T<:Any, A::MPolyRingElem)
    # TODO: test the reconstruction
    # P = [Lq, a, b] encodes sets (x,a(x)/b(x)) where Lq[i](x)=0
    # Compute divisor dA of q that encodes intersection with A(x,y)=0
    if iszero(A)
        return P[1]
    end
    compt = 0
    RS = A |> parent |> symbols
    dA, dAf = [], []
    pprod, p = 1, ZZ(1) << 60
    while compt<12
        #print("$compt,")
        if compt>0
            dAfold, dAold = copy(dAf), copy(dA)
        end
        p = next_prime(p)
        Pp = change_coefficient_ring.(Ref(GF(p)), P)
        Ap = change_coefficient_ring(GF(p), A)
        #@time "Eval"
        Apev = num_biv_rat_mod(Ap, Pp, RS)
        dAp = gcd(Pp[1], Apev)
        dA = lift.(Ref(ZZ), coefficients_of_univariate(dAp))
        if compt>0
            dA = [ crt([d1, d2], [pprod, p], false) for (d1, d2) in zip(dAold, dA) ]
        end
        if length(filter(!is_zero, dA)) <= 1
            #println("Trivial gcd")
            return one(polynomial_ring(QQ, RS)[1])
        end
        pprod = pprod*p
        try
            dAf = [ reconstruct(c, pprod) for c in dA ]
            if compt>0
                dAf != dAfold || break
            end
        catch
            #println("No rational rec..")
        end
        compt += 1
    end
    #println()
    B, = polynomial_ring(QQ, RS[1])
    return int_coeffs(change_ringvar(B(dAf), RS))
end

function param_crit_split(f, g; v=0, detect_app=true)
    # Compute subresultants
    v>0 && println("Compute subresultant sequence")
    @iftime v>0 sr = subresultants(f, derivative(f, 2), 2, list=true)

    if total_degree(sr[1][1]) == 0
        return Dict()
    end

    # TODO: sr is univariate? check what's called
    v>1 && println("Factorization")
    @iftime v>1 sqr = collect(factor(sr[1][1]))

    # Filter out factors with zero multiplicity and sort by multiplicity
    sqr = sort!(filter(t -> t[2] > 0, sqr), by = t -> t[2])
    # Group factors by multiplicity
    sqrmult = unique(getindex.(sqr, 2))
    group_sqr = Dict(m => [r[1] for r in sqr if r[2] == m] for m in sqrmult)

    v>0 && println("Compute crit partition w.r.t to multiplicity")
    # Initalization
    singmult = filter(p->p*(p-1)<=sqrmult[end], 2:sqrmult[end])
    param_crit = Dict(p => [[], -sr[p][end-1], (p-1)*sr[p][end]]  for p in singmult)

    # Critical points : multiplicity 1 in res
    (1 in sqrmult) && push!(param_crit, 1=>[group_sqr[1], -sr[2][end-1], sr[2][end]])
    if length(singmult) == 0
        return filter(p->length(p[2][1])>0, param_crit)
    end

    # Nodes : multiplicity 2 in res
    v>0 && println("Compute apparent singularities")
    if 2 in sqrmult
        if detect_app
        @iftime v>0 begin
            A = derivative(derivative(f,2),2)*derivative(g,1) - derivative(derivative(f,1),2)*derivative(g,2)
            dA = [ (intersect_biv([q, -sr[2][end-1], sr[2][end]], A)) for q in group_sqr[2] ]
            push!(param_crit, 0 => [group_sqr[2]./dA, -sr[2][end-1], sr[2][end]])
            append!(param_crit[2][1], dA)
        end
        else
            push!(param_crit, 0 => [group_sqr[2], -sr[2][end-1], sr[2][end]])
        end
    end
    # Other sing
    filter!(m->!(m in [1,2]), sqrmult)
    lsr = [ try sr[p][end] catch; one(parent(f)) end for p in 2:(singmult[end]+1) ]
    Ld = Vector{Vector{Dict{Int, MPolyRingElem}}}(undef, length(sqrmult))
    v>0 && println("Compute gcd decomposition")
@iftime v>0 begin
    for k in eachindex(sqrmult)
        Ld[k] = Vector{Dict{Int, MPolyRingElem}}(undef,length(group_sqr[sqrmult[k]]))
        for l in eachindex(group_sqr[sqrmult[k]])
            Ld[k][l] = fact_gcd(group_sqr[sqrmult[k]][l], lsr)
        end
    end
    for k in eachindex(Ld)
        for l in eachindex(Ld[k])
            for (i, dji) in Ld[k][l]
                if haskey(param_crit, i+1)
                    push!(param_crit[i+1][1], dji)
                else
                    error("Curve not in generic position")
                end
            end
        end
    end
end
    for i in eachindex(param_crit)
        filter!(p->total_degree(p)>0, param_crit[i][1])
    end
    return filter(p->length(p[2][1])>0, param_crit)
end

param_crit_split(f; v=1) = param_crit_split(f, zero(parent(f)), v=v, detect_app=false)

