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
    
    gb = isempty(values(I.gb)) ? groebner_basis(I) : first(values(I.gb))
    R = parent(first(gb))
    res = [trues(ngens(R))]

    lms = (Nemo.leading_monomial).(gb)
    for lm in lms
        to_del = Int[]
        new_miss = BitVector[]
        for (i, mis) in enumerate(res)
            nz_exps_inds = findall(e -> !iszero(e),
                                   first(Nemo.exponent_vectors(lm)))
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

    max_length = maximum(mis -> length(findall(mis)), res)
    return max_length
end
