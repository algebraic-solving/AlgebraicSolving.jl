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
    
    gb = get(I.gb, 0, groebner_basis(I, complete_reduction = true))
    R = parent(first(gb))
    res = [trues(ngens(R))]

    lead_exps = (_drl_lead_exp).(gb)
    for lexp in lead_exps 
        to_del = Int[]
        new_miss = BitVector[]
        for (i, mis) in enumerate(res)
            nz_exps_inds = findall(e -> !iszero(e), lexp)
            ind_var_inds = findall(mis)
            if issubset(nz_exps_inds, ind_var_inds)
                for j in nz_exps_inds
                    new_mis = copy(mis)
                    new_mis[j] = false
                    push!(new_miss, new_mis)
                end
                push!(to_del, i)
            end
        end
        deleteat!(res, to_del)
        append!(res, new_miss)
        unique!(res)
    end

    length(res) == 0 && return -1
    max_length = maximum(mis -> length(findall(mis)), res)
    return max_length
end

function _drl_exp_vector(u::Vector{Int})
    return [sum(u), -reverse(u)...]
end

function _drl_lead_exp(p::MPolyRingElem)
    exps = collect(Nemo.exponent_vectors(p))
    _, i = findmax((u -> _drl_exp_vector(u)).(exps))
    return exps[i]
end
