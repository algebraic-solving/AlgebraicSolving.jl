@doc Markdown.doc"""
    dimension(I::Ideal{T}) where T <: MPolyRingElem

Compute the Krull dimension of a given polynomial ideal `I`.

**Note**: This requires a Gröbner basis of `I`, which is computed internally if not already known.

# Examples
```jldoctest
julia> using AlgebraicSolving

julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"]);

julia> I = Ideal([x*y,x*z,y*z]);

julia> dimension(I)
1
```
"""
function dimension(I::Ideal{T}) where T <: MPolyRingElem

    !isnothing(I.dim) && return I.dim
    gb = get!(I.gb, 0) do
        groebner_basis(I, complete_reduction = true)
    end
    R = parent(first(gb))

    res = Set([trues(ngens(R))])
    lead_exps = [ _lead_exp_ord(g, :degrevlex) for g in gb if !iszero(g) ]
    for lexp in lead_exps
        nz_exps = (!iszero).(lexp)
        nz_exps_ind = findall(nz_exps)
        next_res = Set{BitVector}()
        for mis in res
            if _all_lesseq(nz_exps, mis)
                @inbounds for j in nz_exps_ind
                    new_mis = copy(mis)
                    new_mis[j] = false
                    push!(next_res, new_mis)
                end
            else
                push!(next_res, mis)
            end
        end
        res = next_res
    end

    I.dim = isempty(res) ? -1 : maximum(sum, res)
    return I.dim
end

function _all_lesseq(a::BitVector, b::BitVector)::Bool
    @inbounds for i in eachindex(a)
        if a[i] && !b[i]
            return false
        end
    end
    return true
end

function _lead_exp_ord(p::MPolyRingElem, order::Symbol)
    @req !iszero(p) "Zero polynomial does not have a leading term"
    R = parent(p)
    internal_ordering(R)==order && return first(exponent_vectors(p))

    A = base_ring(R)
    R1, _ = polynomial_ring(A, symbols(R), internal_ordering=order)
    ctx = MPolyBuildCtx(R1)
    for e in exponent_vectors(p)
        push_term!(ctx, one(A), e)
    end
    return first(exponent_vectors(finish(ctx)))
end