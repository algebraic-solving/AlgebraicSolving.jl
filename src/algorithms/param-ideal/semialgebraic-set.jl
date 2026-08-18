abstract type SemialgebraicSet end

@kwdef struct Empty <: SemialgebraicSet end

@kwdef struct Universe <: SemialgebraicSet end

@kwdef struct BasicSemialgebraicSet <: SemialgebraicSet
    eqs::Vector{QQMPolyRingElem}
    ineqs::Vector{QQMPolyRingElem}
    pos::Vector{QQMPolyRingElem}
    nonneg::Vector{QQMPolyRingElem}
end

@kwdef struct SignEnum <: SemialgebraicSet
    polys::Vector{QQMPolyRingElem}
    signs::Vector{Vector{Int}}
    counts::Vector{Int}
    witnesses::Vector{Vector{QQFieldElem}}
end

@kwdef struct SigEnum <: SemialgebraicSet
    hms::Vector{MatSpaceElem{FracFieldElem{QQMPolyRingElem}}}
    sigs::Vector{Vector{Int}}
    counts::Vector{Int}
    witnesses::Vector{Vector{QQFieldElem}}
end

@kwdef struct MatEnum <: SemialgebraicSet
    hm::MatSpaceElem{FracFieldElem{QQMPolyRingElem}}
    mm::Vector{MatSpaceElem{FracFieldElem{QQMPolyRingElem}}}
end

struct Complement <: SemialgebraicSet
    S::SemialgebraicSet
end

struct Intersect <: SemialgebraicSet
    S::Vector{SemialgebraicSet}
end

struct UnionSet <: SemialgebraicSet
    S::Vector{SemialgebraicSet}
end

# print methods

function Base.show(io::IO, _::Empty)
    print(io, "Empty()")
end

function Base.show(io::IO, _::Universe)
    print(io, "Universe()")
end

function Base.show(io::IO, S::BasicSemialgebraicSet)
    print(io, "BasicSemialgebraicSet(")
    print(io, "eqs=[")
    for (i, eq) in enumerate(S.eqs)
        if i > 1
            print(io, ", ")
        end
        print(io, eq)
    end
    print(io, "], ineqs=[")
    for (i, ineq) in enumerate(S.ineqs)
        if i > 1
            print(io, ", ")
        end
        print(io, ineq)
    end
    print(io, "], pos=[")
    for (i, p) in enumerate(S.pos)
        if i > 1
            print(io, ", ")
        end
        print(io, p)
    end
    print(io, "], nonneg=[")
    for (i, p) in enumerate(S.nonneg)
        if i > 1
            print(io, ", ")
        end
        print(io, p)
    end
    print(io, "])")
end

function Base.show(io::IO, S::SignEnum)
    print(io, "SignEnum(")
    print(io, "polys=[")
    for (i, p) in enumerate(S.polys)
        if i > 1
            print(io, ", ")
        end
        print(io, p)
    end
    print(io, "], signs=")
    print(io, S.signs)
    print(io, ", counts=")
    print(io, S.counts)
    print(io, ", witnesses=")
    print(io, S.witnesses)
    print(io, ")")
end

function Base.show(io::IO, S::SigEnum)
    print(io, "SigEnum(")
    print(io, "hms=[")
    for (i, hm) in enumerate(S.hms)
        if i > 1
            print(io, ", ")
        end
        print(io, hm)
    end
    print(io, "], sigs=")
    print(io, S.sigs)
    print(io, ", counts=")
    print(io, S.counts)
    print(io, ", witnesses=")
    print(io, S.witnesses)
    print(io, ")")
end

function Base.show(io::IO, S::MatEnum)
    print(io, "MatEnum(")
    print(io, "hm=")
    print(io, S.hm)
    print(io, ", mm=[")
    for (i, mm) in enumerate(S.mm)
        if i > 1
            print(io, ", ")
        end
        print(io, mm)
    end
    print(io, "])")
end

function Base.show(io::IO, S::Complement)
    print(io, "Complement(")
    show(io, S.S)
    print(io, ")")
end

function Base.show(io::IO, S::Intersect)
    print(io, "Intersect([")
    for (i, Sᵢ) in enumerate(S.S)
        if i > 1
            print(io, ", ")
        end
        show(io, Sᵢ)
    end
    print(io, "])")
end

function Base.show(io::IO, S::UnionSet)
    print(io, "UnionSet([")
    for (i, Sᵢ) in enumerate(S.S)
        if i > 1
            print(io, ", ")
        end
        show(io, Sᵢ)
    end
    print(io, "])")
end

function Base.show(io::IO, S::SemialgebraicSet)
    if S isa Empty
        show(io, S::Empty)
    elseif S isa Universe
        show(io, S::Universe)
    elseif S isa BasicSemialgebraicSet
        show(io, S::BasicSemialgebraicSet)
    elseif S isa SignEnum
        show(io, S::SignEnum)
    elseif S isa SigEnum
        show(io, S::SigEnum)
    elseif S isa MatEnum
        show(io, S::MatEnum)
    elseif S isa Complement
        show(io, S::Complement)
    elseif S isa Intersect
        show(io, S::Intersect)
    elseif S isa UnionSet
        show(io, S::UnionSet)
    else
        error("Unknown SemialgebraicSet type")
    end
end

# equality methods

function Base.:(==)(_::Empty, _::Empty)
    return true
end

function Base.:(==)(_::Universe, _::Universe)
    return true
end

function Base.:(==)(S1::BasicSemialgebraicSet, S2::BasicSemialgebraicSet)
    return S1.eqs == S2.eqs && S1.ineqs == S2.ineqs && S1.pos == S2.pos && S1.nonneg == S2.nonneg
end

function Base.:(==)(S1::SignEnum, S2::SignEnum)
    return S1.polys == S2.polys && S1.signs == S2.signs && S1.counts == S2.counts
end

function Base.:(==)(S1::SigEnum, S2::SigEnum)
    return S1.hms == S2.hms && S1.sigs == S2.sigs && S1.counts == S2.counts
end

function Base.:(==)(S1::MatEnum, S2::MatEnum)
    return S1.hm == S2.hm && S1.mm == S2.mm
end

function Base.:(==)(S1::Complement, S2::Complement)
    return S1.S == S2.S
end

function Base.:(==)(S1::Intersect, S2::Intersect)
    return S1.S == S2.S
end

function Base.:(==)(S1::UnionSet, S2::UnionSet)
    return S1.S == S2.S
end

function Base.:(==)(S1::SemialgebraicSet, S2::SemialgebraicSet)
    if S1 isa BasicSemialgebraicSet && S2 isa BasicSemialgebraicSet
        return S1 == S2
    elseif S1 isa SignEnum && S2 isa SignEnum
        return S1 == S2
    elseif S1 isa SigEnum && S2 isa SigEnum
        return S1 == S2
    elseif S1 isa MatEnum && S2 isa MatEnum
        return S1 == S2
    elseif S1 isa Complement && S2 isa Complement
        return S1 == S2
    elseif S1 isa Intersect && S2 isa Intersect
        return S1 == S2
    elseif S1 isa UnionSet && S2 isa UnionSet
        return S1 == S2
    else
        return false
    end
end

# evaluation methods

function evaluate(_::Empty, _::Vector{QQFieldElem})::Bool
    return false
end

function evaluate(_::Universe, _::Vector{QQFieldElem})::Bool
    return true
end

function evaluate(S::BasicSemialgebraicSet, vals::Vector{QQFieldElem})::Bool
    for eq in S.eqs
        if evaluate(eq, vals) != 0
            return false
        end
    end
    for ineq in S.ineqs
        if evaluate(ineq, vals) == 0
            return false
        end
    end
    for p in S.pos
        if evaluate(p, vals) <= 0
            return false
        end
    end
    for p in S.nonneg
        if evaluate(p, vals) < 0
            return false
        end
    end
    return true
end

function evaluate(S::SignEnum, vals::Vector{QQFieldElem})::Bool
    signs = [sign(evaluate(p, vals)) for p in S.polys]
    i = findfirst(==(signs), S.signs)
    return isnothing(i) ? false : S.counts[i] > 0
end

function evaluate(S::SigEnum, vals::Vector{QQFieldElem})::Bool
    hm_specs = map(hm -> map(a -> evaluate(numerator(a), vals) // evaluate(denominator(a), vals), hm), S.hms)
    sig = map(hm -> signature(hm), hm_specs)
    i = findfirst(==(sig), S.sigs)
    return isnothing(i) ? false : S.counts[i] > 0
end

function evaluate(S::MatEnum, vals::Vector{QQFieldElem})::Bool
    function _sum_of_signatures(_S::MatEnum, _vals::Vector{QQFieldElem})::Int
        if isempty(_S.mm)
            return signature(map(a -> evaluate(numerator(a), _vals) // evaluate(denominator(a), _vals), _S.hm))
        end
        total = _sum_of_signatures(MatEnum(hm=_S.hm, mm=_S.mm[2:end]), _vals)
        total += _sum_of_signatures(MatEnum(hm=_S.hm * _S.mm[1], mm=_S.mm[2:end]), _vals)
        return total
    end
    return _sum_of_signatures(S, vals) > 0
end

function evaluate(S::Complement, vals::Vector{QQFieldElem})::Bool
    return !evaluate(S.S, vals)
end

function evaluate(S::Intersect, vals::Vector{QQFieldElem})::Bool
    for Sᵢ in S.S
        if !evaluate(Sᵢ, vals)
            return false
        end
    end
    return true
end

function evaluate(S::UnionSet, vals::Vector{QQFieldElem})::Bool
    for Sᵢ in S.S
        if evaluate(Sᵢ, vals)
            return true
        end
    end
    return false
end

function evaluate(S::SemialgebraicSet, vals::Vector{QQFieldElem})::Bool
    if S isa Empty
        return evaluate(S::Empty, vals)
    elseif S isa Universe
        return evaluate(S::Universe, vals)
    elseif S isa BasicSemialgebraicSet
        return evaluate(S::BasicSemialgebraicSet, vals)
    elseif S isa SignEnum
        return evaluate(S::SignEnum, vals)
    elseif S isa SigEnum
        return evaluate(S::SigEnum, vals)
    elseif S isa MatEnum
        return evaluate(S::MatEnum, vals)
    elseif S isa Complement
        return evaluate(S::Complement, vals)
    elseif S isa Intersect
        return evaluate(S::Intersect, vals)
    elseif S isa UnionSet
        return evaluate(S::UnionSet, vals)
    else
        error("Unknown SemialgebraicSet type")
    end
end
