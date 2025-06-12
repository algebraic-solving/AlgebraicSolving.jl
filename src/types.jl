export Ideal

mutable struct RationalParametrization
    vars::Vector{Symbol}
    cfs_lf::Vector{ZZRingElem}
    elim::QQPolyRingElem
    denom::QQPolyRingElem
    param::Vector{QQPolyRingElem}

    function RationalParametrization(
            vars::Vector{Symbol},
            cfs_lf::Vector{ZZRingElem},
            elim::QQPolyRingElem,
            denom::QQPolyRingElem,
            param::Vector{QQPolyRingElem}
        )
        rp = new()
        rp.vars   = vars
        rp.cfs_lf = cfs_lf
        rp.elim   = elim
        rp.denom  = denom
        rp.param  = param

        return rp
    end
end

mutable struct RationalCurveParametrization
    vars::Vector{Symbol}
    cfs_lfs::Vector{Vector{ZZRingElem}}
    elim::QQMPolyRingElem
    denom::QQMPolyRingElem
    param::Vector{QQMPolyRingElem}

    function RationalCurveParametrization(
            vars::Vector{Symbol},
            cfs_lfs::Vector{Vector{ZZRingElem}},
            elim::QQMPolyRingElem,
            denom::QQMPolyRingElem,
            param::Vector{QQMPolyRingElem}
        )
        rp = new()
        rp.vars   = vars
        rp.cfs_lfs = cfs_lfs
        rp.elim   = elim
        rp.denom  = denom
        rp.param  = param

        return rp
    end
end

mutable struct Ideal{T <: MPolyRingElem}
    gens::Vector{T}
    dim::Union{Int, Nothing}
    deg::Union{Int, Nothing}
    gb::Dict{Int, Vector{T}}
    inter_sols::Vector{Vector{Vector{QQFieldElem}}}
    real_sols::Vector{Vector{QQFieldElem}}
    rat_sols::Vector{Vector{QQFieldElem}}
    rat_param::RationalParametrization

    function Ideal(F::Vector{T}) where {T <: MPolyRingElem}
        I = new{T}()
        I.gens = F
        I.gb   = Dict()
        I.dim  = nothing
        I.deg  = nothing
        return I
    end
end

Base.parent(I::Ideal) = Nemo.parent(I.gens[1])

Base.show(io::IO, I::Ideal) = print(io, I.gens)

Base.getindex(I::Ideal, idx::Int) = I.gens[idx]
