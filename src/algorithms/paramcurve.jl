export compute_param, param_newvars


# Return the polynomials in F, but injected in the polynomial ring with newvarias_S as new variables
function change_ringvar(F::Vector{P}, newvarias_S::Vector{Symbol}) where {P <: MPolyRingElem}
    R = parent(first(F))
    # Locate variables of R in newvarias
    to_varias = Vector{Int}(undef,0)
    for v in newvarias_S
        ind = findfirst(x->x==v, R.S)
        push!(to_varias, typeof(ind)==Nothing ? length(R.S)+1 : ind)
    end

    ind_novarias = setdiff(eachindex(R.S), to_varias)
    newR, newvarias = polynomial_ring(base_ring(R), newvarias_S)

    res = typeof(first(F))[]
    ctx = MPolyBuildCtx(newR)

    for f in F
        for (e, c) in zip(exponent_vectors(f), coefficients(f))
            @assert(all([ e[i]==0 for i in ind_novarias ]), "Occurence of old variable.s found!")
            push!(e, 0)
            push_term!(ctx, c, [e[i] for i in to_varias ])
        end
        push!(res, finish(ctx))
    end

    return res
end

function deg_Alg(F, dim)
    if dim <= 0
        r = rational_parametrization(Ideal(F))
    else
        varias = gens(parent(first(F)))
        planes = [ sum([ rand(ZZ(-100):ZZ(100))*v for v in vcat(varias,1) ])  for _ in 1:dim]
        r = rational_parametrization(Ideal(vcat(F,planes)))
    end
    return degree(r.elim)
end

function compute_param(I::Ideal{P} where P<:MPolyRingElem; use_lfs = false, 
    cfs_lfs::Vector{Vector{ZZRingElem}} = Vector{Vector{ZZRingElem}}())
    R = parent(I)
    varias, N = gens(R), nvars(R)
    F = I.gens
    DEG = deg_Alg(F,1)

    nb_lfs = length(cfs_lfs)
    @assert(nb_lfs in [0,2], "Provide either none or two linear forms")
    @assert(all(length(cfs_lf) == N+2 for cfs_lf in cfs_lfs), "Linear form(s) of wrong size")
    if use_lfs || nb_lfs>0
        # Add new generic variables at the end if necessary
        newvarias_S = vcat(R.S, [:B,:A])
        # We change the polynomial ring and add linear form(s) if necessary
        R, varias = polynomial_ring(base_ring(R), newvarias_S)
        N = length(varias)
        # Add linear forms to keep dimension
        F = change_ringvar(F, newvarias_S)
        if nb_lfs == 0
            cfs_lfs = [rand(ZZ.(setdiff(-100:100,0)), N) for _ in 1:2]
        end
        append!(F, [ transpose(cfs_lf)*varias  for cfs_lf in cfs_lfs ])
    end

    @assert(deg_Alg(vcat(F, varias[N-1]-rand(ZZ(-100):ZZ(100))), 0) == DEG,
         "The curve is not in generic position")

    # Compute DEG+1 evaluations of the param (whose deg is bounded by DEG)
    PARAM  = Vector{Vector{AlgebraicSolving.QQPolyRingElem}}(undef,DEG+2)
    _values = Vector{QQFieldElem}(undef,DEG+2)
    i = 1
    free_ind = collect(1:DEG+2)
    used_ind = zeros(Bool, DEG+2)
    while length(free_ind) > 0
        if i > 2*(DEG+2)
            error("Too many bad specializations")
        end
        # Evaluation of the generators
        LFeval = Vector{AlgebraicSolving.Ideal}(undef, length(free_ind))
        Threads.@threads for j in 1:length(free_ind)
            LFeval[i+j-1] = Ideal(change_ringvar(evaluate.(F, Ref([N-1]), Ref([QQ(i+j-1)])), [R.S[1:N-2]; R.S[N]]))
        end
        # Compute parametrization of each evaluation
        Lr = Vector{AlgebraicSolving.RationalParametrization}(undef, length(free_ind))
        for j in 1:length(free_ind)
            Lr[i+j-1] = rational_parametrization(LFeval[i+j-1], nr_thrds=Threads.nthreads())
        end
        for j in 1:length(free_ind)
            # For lifting: the same variable must be chosen for the param
            if  Lr[j].vars == [R.S[1:N-2]; R.S[N]]
                lc = leading_coefficient(Lr[j].elim)
                rr = [ p/lc for p in vcat(Lr[j].elim, Lr[j].denom, Lr[j].param) ]
                PARAM[j] = rr
                _values[j] = QQ(i+j-1)
                used_ind[j] = true
            else
                println("bad specialization: ", i+j-1)
            end
        end

        i += length(free_ind)
        free_ind = [ free_ind[j] for j in eachindex(free_ind) if !used_ind[j] ]
        used_ind = zeros(Bool, length(free_ind))
    end

    # Interpolate each coefficient of each poly in the param
    POLY_PARAM = Vector{QQMPolyRingElem}(undef,N)
    T, (x,y) = polynomial_ring(QQ, [:x,:y])
    A, = polynomial_ring(QQ, :u)
    Threads.@threads for count in 1:N
        COEFFS = Vector{QQPolyRingElem}(undef, DEG+1)
        Threads.@threads for deg in 0:DEG
            _evals = [coeff(PARAM[i][count], deg) for i in 1:length(PARAM)]
            COEFFS[deg+1] = interpolate(A, _values, _evals)
        end

        C = [ collect(coefficients(c)) for c in COEFFS ]
        POL_term = [C[i][j]*y^(i-1)*x^(j-1) for i in eachindex(C) for j in eachindex(C[i])]
        POL = length(POL_term) > 0 ? sum(POL_term) : T(0)

        POLY_PARAM[count] = POL
    end

    return [R.S, cfs_lfs, POLY_PARAM[1], POLY_PARAM[2], POLY_PARAM[3:end]]
end
