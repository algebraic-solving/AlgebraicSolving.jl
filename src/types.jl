export Ideal

mutable struct RationalParametrization
    vars::Vector{Symbol}
    lf_cfs::Vector{ZZRingElem}
    elim::QQPolyRingElem
    denom::QQPolyRingElem
    param::Vector{PolyRingElem}

    function RationalParametrization(
            vars::Vector{Symbol},
            lf_cfs::Vector{ZZRingElem},
            elim::QQPolyRingElem,
            denom::QQPolyRingElem,
            param::Vector{PolyRingElem}
        )
        rp = new()
        rp.vars   = vars
        rp.lf_cfs = lf_cfs
        rp.elim   = elim
        rp.denom  = denom
        rp.param  = param

        return rp
    end
end

mutable struct Ideal{T <: MPolyRingElem}
    gens::Vector{T}
    dim::Int
    gb::Dict{Int, Vector{T}}
    inter_sols::Vector{Vector{Vector{QQFieldElem}}}
    real_sols::Vector{Vector{QQFieldElem}}
    rat_sols::Vector{Vector{QQFieldElem}}
    rat_param::RationalParametrization

    function Ideal(F::Vector{T}) where {T <: MPolyRingElem}
        I = new{T}()
        I.gens = F
        I.dim  = -1
        I.gb   = Dict()

        return I
    end
end

Base.parent(I::Ideal) = Nemo.parent(I.gens[1])

Base.show(io::IO, I::Ideal) = print(io, I.gens)

Base.getindex(I::Ideal, idx::Int) = I.gens[idx]
