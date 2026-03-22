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

function change_ringvar(F::Vector{P}, newvarias_S::Vector{Symbol}) where {P <: PolyRingElem}
    R = parent(first(F))
    to_varias = [ v==R.S ? 1 : 2 for v in newvarias_S]
    newR, = polynomial_ring(base_ring(R), newvarias_S)

    res = typeof(zero(newR))[]
    ctx = MPolyBuildCtx(newR)

    LcF = [ filter!(t->t[2]!=0, collect(enumerate(coefficients(f)))) for f in F ]
    for f in LcF
        for (ef, c) in f
            e = [ef-1, 0]
            push_term!(ctx, c, [e[i] for i in to_varias ])
        end
        push!(res, finish(ctx))
    end

    return res
end

function change_ringvar(f::Union{MPolyRingElem, PolyRingElem}, newvarias_S::Vector{Symbol})
    return first(change_ringvar([f], newvarias_S))
end

function change_ringvar(F::Vector{P}, ind_newvarias::Union{Vector{I}, UnitRange{I}}) where {P <: MPolyRingElem, I <: Int64}
    R = parent(first(F))
    return change_ringvar(F, [R.S[i] for i in ind_newvarias])
end

function change_ringvar(f::MPolyRingElem, ind_newvarias::Union{Vector{I}, UnitRange{I}}) where {I <: Int64}
    R = parent(f)
    return first(change_ringvar([f], [R.S[i] for i in ind_newvarias]))
end

# Return the polynomials in F, but injected in the polynomial ring with the variables occuring in F
function change_ringvar(F::Vector{P}) where {P <: MPolyRingElem}
    union_varias = Set{Symbol}()
    for f in F
        union!(union_varias, map(Symbol, vars(f)) )
    end
    return change_ringvar(F, collect(union_varias))
end

function change_ringvar(f::MPolyRingElem)
    return first(change_ringvar([f]))
end

# Return the polynomials in F, but injected in the polynomial ring with newvarias_S as new variables
function change_ringvar_mod(F::Vector{P}, newvarias_S::Vector{Symbol}, oldvarias_S::Vector{Symbol}) where {P <: MPolyRingElem}
    R = parent(first(F))
    # Locate variables of R in newvarias
    to_varias = Vector{Int}(undef,0)
    for v in newvarias_S
        ind = findfirst(x->x==v, oldvarias_S)
        push!(to_varias, typeof(ind)==Nothing ? length(oldvarias_S)+1 : ind)
    end

    ind_novarias = setdiff(eachindex(oldvarias_S), to_varias)
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

function change_ringvar_mod(F::Vector{P}, newvarias_S::Vector{Symbol}, oldvarias_S::Vector{Symbol}) where {P <: PolyRingElem}
    R = parent(first(F))
    A, (x,) = polynomial_ring(base_ring(R), oldvarias_S)
    return change_ringvar_mod([ evaluate(f, x) for f in F ], newvarias_S, oldvarias_S)
end

function change_ringvar_mod(f::Union{MPolyRingElem, PolyRingElem}, newvarias_S::Vector{Symbol}, oldvarias_S::Vector{Symbol})
    return first(change_ringvar_mod([f], newvarias_S, oldvarias_S))
end

#function change_ringvar_mod(f::Union{MPolyRingElem, PolyRingElem}, newvarias_S::Vector{Symbol})
#    return first(change_ringvar_mod([f], newvarias_S, oldvarias_S))
#end

function MPolyBuild(F::Vector{Vector{P}}, newvarias_S::Vector{Symbol}, idx::Int) where {P <: RingElem}
    A = parent(first(first(F)))
    R, = polynomial_ring(A, newvarias_S)
    to_varias = [ i==idx ? 1 : 2 for i in eachindex(newvarias_S) ]

    res = typeof(zero(R))[]
    ctx = MPolyBuildCtx(R)

    LcF = [ filter!(t->t[2]!=0, collect(enumerate(f))) for f in F ]
    for f in LcF
        for (ef, c) in f
            e = [ef-1, 0]
            push_term!(ctx, c, [e[i] for i in to_varias ])
        end
        push!(res, finish(ctx))
    end

    return res
end

function MPolyBuild(f::Vector{P}, newvarias_S::Vector{Symbol}, idx::Int) where {P <: RingElem}
    return first(MPolyBuild([f], newvarias_S, idx))
end