export compute_param, add_genvars

function compute_param(
        I::Ideal{P} where P<:MPolyRingElem;                         # input generators
        use_lfs::Bool = false,                                      # add generic variables
        cfs_lfs::Vector{Vector{ZZRingElem}} = Vector{ZZRingElem}[]  # coeffs of linear forms
    )
    Itest = Ideal(change_base_ring.(Ref(GF(65521)), I.gens))
    @assert(I.dim==1 || I.dim<0 && dimension(Itest)==1, "I must be one-dimensional dimension")
    DEG = hilbert_degree(Itest)
    # If generic variables must be added
    F, cfs_lfs = use_lfs ? add_genvars(I.gens, max(2, length(cfs_lfs)), cfs_lfs) :
             !isempty(cfs_lfs) ? add_genvars(I.gens, length(cfs_lfs), cfs_lfs) :
             (I.gens, Vector{ZZRingElem}[])
    R = parent(first(F))
    N = nvars(R)
    let # local code block test
        local INEW = change_base_ring.(Ref(GF(65521)), vcat(F, gens(R)[N-1]-rand(ZZ,-100:100))) |> Ideal
        @assert(dimension(INEW)==0 && hilbert_degree(INEW) == DEG, "The curve is not in generic position")
    end

    # Compute DEG+2 evaluations of x in the param (whose total deg is bounded by DEG)
    PARAM  = Vector{Vector{AlgebraicSolving.QQPolyRingElem}}(undef,DEG+2)
    _values = Vector{QQFieldElem}(undef,DEG+2)
    i = 1
    free_ind = collect(1:DEG+2)
    used_ind = zeros(Bool, DEG+2)
    while length(free_ind) > 0
        if i > 2*(DEG+2)
            error("Too many bad specializations: permute variables or use_lfs=true")
        end
        # Evaluation of the generators
        LFeval = Vector{AlgebraicSolving.Ideal}(undef, length(free_ind))
        Threads.@threads for j in 1:length(free_ind)
            LFeval[j] = Ideal(change_ringvar(evaluate.(F, Ref([N-1]), Ref([QQ(i+j-1)])), [R.S[1:N-2]; R.S[N]]))
        end
        # Compute parametrization of each evaluation
        Lr = Vector{AlgebraicSolving.RationalParametrization}(undef, length(free_ind))
        for j in 1:length(free_ind)
            Lr[j] = rational_parametrization(LFeval[j], nr_thrds=Threads.nthreads())
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
    T = polynomial_ring(QQ, [:x,:y])[1]
    A = polynomial_ring(QQ)[1]

    POLY_PARAM = Vector{QQMPolyRingElem}(undef,N)
    for count in 1:N
        COEFFS = Vector{QQPolyRingElem}(undef, DEG+1)
        Threads.@threads for deg in 0:DEG
            _evals = [coeff(PARAM[i][count], deg) for i in 1:length(PARAM)]
            COEFFS[deg+1] = interpolate(A, _values, _evals)
        end
        ctx = MPolyBuildCtx(T)
        for (i, c) in enumerate(COEFFS)
            for (j, coeff) in enumerate(coefficients(c))
                !iszero(coeff) && push_term!(ctx, coeff, [j-1, i-1])
            end
        end
        POLY_PARAM[count] = finish(ctx)
    end

    # Output: [vars, linear forms, elim, denom, [nums_param]]
    return RationalCurveParametrization(R.S, cfs_lfs, POLY_PARAM[1],
                                        POLY_PARAM[2], POLY_PARAM[3:end])
end


# Return the polynomials in F, but injected in the polynomial ring with newvarias_S as new variables
function change_ringvar(
        F::Vector{P} where P <: MPolyRingElem,  # list of polynomials
        newvarias_S::Vector{Symbol}             # new variable symbols
        )
    R = parent(first(F))
    # Locate variables of R in newvarias
    to_varias = Vector{Int}(undef,0)
    for v in newvarias_S
        ind = findfirst(x->x==v, R.S)
        push!(to_varias, typeof(ind)==Nothing ? length(R.S)+1 : ind)
    end

    ind_novarias = setdiff(eachindex(R.S), to_varias)
    newR, = polynomial_ring(base_ring(R), newvarias_S)

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

function add_genvars(
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
newS = vcat(R.S, Symbol.(["_Z$i" for i in ngenvars:-1:1]))
R_ext, all_vars = polynomial_ring(base_ring(R), newS)
Fnew = change_ringvar(F, newS)

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

