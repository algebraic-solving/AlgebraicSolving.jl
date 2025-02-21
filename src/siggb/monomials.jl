#-- Monomial arithmetic --#

function mul(a::Monomial, b::Monomial)
    return Monomial(a.deg + b.deg, a.exps + b.exps)
end

function mul!(buf::MVector{N, Exp}, a::Monomial{N}, b::Monomial{N}) where N
    @inbounds for i in 1:N
        buf[i] = a.exps[i] + b.exps[i]
    end
end

function divide(a::Monomial, b::Monomial)
    return Monomial(a.deg - b.deg, a.exps - b.exps)
end

# lcm(a, b) / a
function lcm_div(a::Monomial{N}, b::Monomial{N}) where N
    e1 = a.exps
    e2 = b.exps
    @inbounds exps = SVector{N, Exp}([e1[i] >= e2[i] ? zero(Exp) : e2[i] - e1[i] for i in 1:N])
    return Monomial(sum(exps), exps)
end

# false if the monomial corresponding to a does not divide
# the monomial corresponding to b
@inline function divch(a::DivMask, b::DivMask)
    return iszero(a & (~b))
end

@inline function divch(a::Monomial{N}, b::Monomial{N}) where N
    if a.deg <= b.deg
        @inbounds for i in 1:N
            a.exps[i] > b.exps[i] && return false
        end
        return true
    end
    return false
end

@inline function divch(a::Monomial, b::Monomial,
                       am::DivMask, bm::DivMask)

    return divch(am, bm) && divch(a, b)
end

# TODO: unroll this loop
@inline function divch!(buf::MVector{N, Exp},
                        a::Monomial{N},
                        b::Monomial{N}) where N

    @inbounds for i in 1:N
        buf[i] = a.exps[i] - b.exps[i]
        if buf[i] < 0
            return false
        end
    end
    return true
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

function cmp_ind(ind1::SigIndex, ind2::SigIndex,
                 ind_order::IndOrder)

    ord_ind1 = ind_order.ord[ind1]
    ord_ind2 = ind_order.ord[ind2]
    return ord_ind1 <= ord_ind2
end

function cmp_ind_str(ind1::SigIndex, ind2::SigIndex,
                     ind_order::IndOrder)

    ord_ind1 = ind_order.ord[ind1]
    ord_ind2 = ind_order.ord[ind2]
    return ord_ind1 < ord_ind2
end

function are_incompat(ind1::SigIndex, ind2::SigIndex,
                      ind_order::IndOrder)

    ind1 == ind2 && return false
    k = ind1 < ind2 ? (ind1, ind2) : (ind2, ind1)
    return get(ind_order.incompat, k, false)
end

function lt_pot(a::Sig, b::Sig,
                ind_order::IndOrder)
    if index(a) == index(b)
        return lt_drl(a[2], b[2])
    else
        ord_ind_a = ind_order.ord[index(a)]
        ord_ind_b = ind_order.ord[index(b)]
        return ord_ind_a < ord_ind_b
    end
end

#---------------------#

#-- Masking and hashing --#

# from groebner.jl
function divmask(e::Monomial{N},
                 divmap,
                 ndivbits) where N

    ctr = one(DivMask)
    res = zero(DivMask)
    o = one(DivMask)
    lb = N > DivMaskSize ? N - DivMaskSize : 1
    for i in N:-1:1
        for j in 1:ndivbits
            @inbounds if e.exps[i] >= divmap[ctr]
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
        hash = :($hash + $(rand(rng, MonHash))*(m[$i]))
    end
    return :($hash % MonHash)
end

Base.hash(a::Monomial{N}) where N = makehash(Val(N), a.exps)

function leading_monomial(basis::Basis,
                          basis_ht::MonomialHashtable,
                          i)

    return basis_ht.exponents[first(basis.monomials[i])]
end

#---------------------#

#-- For Polynomials --#

function is_one(pol::Polynomial, ht::MonomialHashtable)
    return length(pol[2]) == 1 && all(iszero, ht.exponents[pol[2][1]].exps)
end

function sort_poly!(pol::Polynomial; kwargs...)
    s = sortperm(pol[2]; kwargs...)
    permute!(pol[1], s)
    permute!(pol[2], s)
    return
end

function normalize_cfs!(cfs::Vector{Coeff},
                        char::Val{Char}) where Char

    isone(first(cfs)) && return
    inver = inv(first(cfs), char)
    @inbounds for i in eachindex(cfs)
        if isone(i)
            cfs[i] = one(Coeff)
            continue
        end
        cfs[i] = mul(inver, cfs[i], char)
    end
end

# assumes mons are sorted ascendingly by hash index
function add_pols(coeffs1::Vector{Coeff},
                  mons1::Vector{MonIdx},
                  coeffs2::Vector{Coeff},
                  mons2::Vector{MonIdx},
                  vch::Val{Char},
                  scalar1::Coeff=one(Coeff),
                  scalar2::Coeff=one(Coeff)) where {Char}

    l1 = length(mons1)
    l2 = length(mons2)
    mons_res = Vector{MonIdx}(undef, l1 + l2)
    coeffs_res = Vector{Coeff}(undef, l1 + l2)
    
    ind1 = 1
    ind2 = 1
    new_l = 0 
    @inbounds while ind1 <= l1 && ind2 <= l2
        new_l += 1
        m1 = mons1[ind1]
        m2 = mons2[ind2]
        if m1 == m2
            mons_res[new_l] = m1
            coeffs_res[new_l] = add(mul(scalar1, coeffs1[ind1], vch),
                                    mul(scalar2, coeffs2[ind2], vch),
                                    vch)
            ind1 += 1
            ind2 += 1
        elseif m1 < m2
            mons_res[new_l] = m1
            coeffs_res[new_l] = mul(scalar1, coeffs1[ind1], vch)
            ind1 += 1
        else
            mons_res[new_l] = m2
            coeffs_res[new_l] = mul(scalar2, coeffs2[ind2], vch)
            ind2 += 1
        end
    end

    while ind1 <= l1
        new_l += 1
        mons_res[new_l] = mons1[ind1]
        coeffs_res[new_l] = mul(scalar1, coeffs1[ind1], vch)
        ind1 += 1
    end

    while ind2 <= l2 
        new_l += 1
        mons_res[new_l] = mons2[ind2]
        coeffs_res[new_l] = mul(scalar2, coeffs2[ind2], vch)
        ind2 += 1
    end

    resize!(mons_res, new_l)
    resize!(coeffs_res, new_l)

    zero_cfs_inds = findall(c -> iszero(c), coeffs_res)
    deleteat!(mons_res, zero_cfs_inds)
    deleteat!(coeffs_res, zero_cfs_inds)

    return coeffs_res, mons_res
end

function mult_pols(exps1::Vector{Monomial{N}},
                   exps2::Vector{Monomial{N}},
                   cfs1::Vector{Coeff},
                   cfs2::Vector{Coeff},
                   char::Val{Char}) where {N, Char}

    R, vrs = polynomial_ring(GF(Int(Char)), ["x$i" for i in 1:N],
                             ordering = :degrevlex)
    p1 = convert_to_pol(R, exps1, cfs1)
    p2 = convert_to_pol(R, exps2, cfs2)
    p = p1*p2

    lp = length(p)
    exps = exponent_vectors(p)
    cfs = coefficients(p)
    
    res_exps = Vector{Monomial{N}}(undef, lp)
    res_cfs = Vector{Coeff}(undef, lp)
    @inbounds for (i, (cf, evec)) in enumerate(zip(cfs, exps)) 
        m = monomial(SVector{N}((Exp).(evec)))
        cff = cf.data
        res_exps[i] = m
        res_cfs[i] = cff
    end

    return res_exps, res_cfs
end

function mul_by_mon(mons::Vector{M},
                    mon::M) where {M <: Monomial}

    mons_res = Vector{M}(undef, length(mons))
    @inbounds for i in 1:length(mons)
        mons_res[i] = mul(mon, mons[i])
    end
    return mons_res
end

function mul_by_coeff(coeffs::Vector{Coeff},
                      c::Coeff,
                      vchar::Val{Char}) where Char 

    coeffs_res = Vector{Coeff}(undef, length(coeffs))
    @inbounds for i in 1:length(coeffs)
        coeffs_res[i] = mul(c, coeffs[i], vchar)
    end
    return coeffs_res
end

function mul_by_coeff!(coeffs::Vector{Coeff},
                       c::Coeff,
                       vchar::Val{Char}) where Char 

    @inbounds for i in 1:length(coeffs)
        coeffs[i] = mul(c, coeffs[i], vchar)
    end
end

#-------------------------#

# for readibility
index(a::Sig) = @inbounds a[1]
monomial(a::Sig) = @inbounds a[2]
index(a::MaskSig) = @inbounds a[1]
mask(a::MaskSig) = @inbounds a[2]

function one_monomial(::Type{Monomial{N}}) where N
    return Monomial{N}(zero(Exp), SVector{N}(zeros(Exp, N)))
end

function Base.show(io::IO, s::Sig)
    show(io, (Int(index(s)), Vector{Int}(monomial(s).exps)))
end
