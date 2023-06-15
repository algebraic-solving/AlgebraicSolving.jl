using StaticArrays

#-- Monomial arithmetic --#

function mul(a::Monomial, b::Monomial)
    return Monomial(a.deg + b.deg, a.exps + b.exps)
end

function divide(a::Monomial, b::Monomial)
    return Monomial(a.deg - b.deg, a.exps - b.exps)
end

# lcm(a, b) / a
function lcm_div(a::Monomial{N}, b::Monomial{N}) where N
    e1 = a.exps
    e2 = b.exps
    @inbounds exps = [e1[i] >= e2[i] ? zero(Exp) : e2[i] - e1[i] for i in 1:N]
    return Monomial(sum(exps), exps)
end

# false if the monomial corresponding to a does not divide
# the monomial corresponding to b
@inline function div(a::DivMask, b::DivMask)
    return iszero(a & (~b))
end

@inline function div(a::Monomial, b::Monomial)
    return a.deg <= b.deg && (a.exps .<= b.exps)
end

# TODO: inlining should make memory access nice
@inline function div(a::Monomial, b::Monomial,
                     am::DivMask, bm::DivMask)

    return div(am, bm) && div(a, b)
end

#-- Monomial Orders --#

# DRL comparison
@generated function lt_drl(a::Monomial{N}, b::Monomial{N}) where N
    quote
        a.deg < b.deg && return true
        a.deg > b.deg && return false
        
        ae = a.exps
        be = b.exps
        $([:(ae[$i] < be[$i] && return false ; ae[$i] > be[$i] && return true) for i in N:-1:1]...)

        return false
    end
end

function lt_pot(a::Sig, b::Sig)
    if index(a) == index(b)
        return lt_drl(a[2], b[2])
    else
        return index(a) < index(b)
    end
end

#---------------------#

#-- Masking and hashing --#

# from groebner.jl
# TODO: make this generated?
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

@generated function makehash(::Val{N}, m) where {N}
    rng = MersenneTwister(18121987)
    hash = :( $(UInt64(0)) )
    for i in 1:N
        hash = :($hash + $(rand(rng, UInt))*UInt(m[$i]))
    end
    return hash
end

Base.hash(a::Monomial{N}) = makehash(Val(N), a.exps)

#-------------------------#

# for readibility
index(a::Sig) = a[1]
index(a::MaskSig) = a[1]
monomial(a::Sig) = a[2]
