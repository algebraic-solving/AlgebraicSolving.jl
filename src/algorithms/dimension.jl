@doc Markdown.doc"""
    dimension(I::Ideal{T}) where T <: MPolyRingElem

Compute the Krull dimension of a given polynomial ideal `I`.

**Note**: This requires a Gröbner basis of `I`, which is computed internally if not alraedy known.

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

    gb = get!(I.gb, 0) do
        groebner_basis(I, complete_reduction = true)
    end
    R = parent(first(gb))

    res = Set([trues(ngens(R))])
    lead_exps = (_drl_lead_exp).(gb)
    for lexp in lead_exps
        nz_exps = (!iszero).(lexp)
        nz_exps_ind = findall(nz_exps)
        next_res = Set{BitVector}()
        for mis in res
            if nz_exps <= mis
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

function _drl_exp_vector(u::Vector{Int})
    return [sum(u), -reverse(u)...]
end

function _drl_lead_exp(p::MPolyRingElem)
    exps = collect(Nemo.exponent_vectors(p))
    _, i = findmax(_drl_exp_vector, exps)
    return exps[i]
end