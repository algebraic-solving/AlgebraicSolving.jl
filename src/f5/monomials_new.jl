using StaticArrays
const Exp = Int16
const DivMask = UInt32
# first entry is degree
const Monomial{N} = SVector{N, Exp}

struct Monomial{N}
    deg::Exp
    exps::SVector{N, Exp}
end

# DRL comparison
@generated function lt_drl(a::Monomial{N}, b::Monomial{N}) where N
    quote
        a[1] < b[1] && return true
        a[1] > b[1] && return false
        
        $([:(a[$i] < b[$i] && return false ; a[$i] > b[$i] && return true) for i in N:-1:2]...)

        return false
    end
end

# from groebner.jl
function divmask(e::Monomial{N},
                 divmap,
                 ndivbits) where N

    ctr = one(DivMask)
    res = zero(DivMask)
    o = one(DivMask)
    for i in 1:N
        for j in 1:ndivbits
            @inbounds if e.exps[i + 1] >= divmap[ctr]
                res |= o << (ctr - 1)
            end
            ctr += o
        end
    end

    res
end

# for readibility

function mul(a::Monomial, b::Monomial)
    return Monomial(a.deg + b.deg, a.exps + b.exps)
end

function divide(a::Monomial, b::Monomial)
    return Monomial(a.deg - b.deg, a.exps - b.exps)
end

function lcm_div(a::Monomial{N}, b::Monomial{N}) where N
    e1 = a.exps
    e2 = b.exps
    @inbounds exps = [e1[i] >= e2[i] ? zero(Exp) : e2[i] - e1[i] for i in 1:N]
    return Monomial(sum(exps), exps)
end

function mask_div(a::DivMask, b::DivMask)
    return iszero(a & (~b))
end

function div(a::Monomial, b::Monomial)
    return a.deg <= b.deg && (a.exps .<= b.exps)
end
