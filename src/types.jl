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
    rat_param::Union{RationalParametrization, RationalCurveParametrization}

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

Base.getindex(I::Ideal, idx::Union{Int, UnitRange}) = I.gens[idx]

Base.lastindex(I::Ideal) = lastindex(I.gens)

mutable struct RMnode
    polar_eqs::Vector{QQMPolyRingElem}
    base_pt::Vector{QQFieldElem}
    children::Vector{RMnode}
end

mutable struct Roadmap
    initial_ideal::Ideal{QQMPolyRingElem}
    root::RMnode
end

function _collect_roadmap(RMn::RMnode, F)
    data = [F(RMn)]
    for child in RMn.children
        append!(data, _collect_roadmap(child, F))
    end
    return data
end

function all_eqs(RM::Roadmap)
    func(s) = fbr(vcat(RM.initial_ideal.gens, s.polar_eqs), s.base_pt)
    return _collect_roadmap(RM.root, func)
end

function all_base_pts(RM::Roadmap)
    return _collect_roadmap(RM.root, s->s.base_pt)
end

function nb_nodes(RM::Roadmap)
    return length(_collect_roadmap(RM.root, s -> true))
end

Base.show(io::IO, RM::Roadmap) = print(io, all_base_pts(RM))
Base.getindex(RM::Roadmap, idx::Union{Int, UnitRange}) = all_eqs(RM)[idx]
Base.lastindex(RM::Roadmap) = nb_nodes(RM)
Base.length(RM::Roadmap) = nb_nodes(RM)