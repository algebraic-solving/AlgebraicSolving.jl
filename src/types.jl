export Ideal, ParametricIdeal

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

mutable struct CurveRationalParametrization
    vars::Vector{Symbol}
    cfs_lfs::Vector{Vector{ZZRingElem}}
    elim::QQMPolyRingElem
    denom::QQMPolyRingElem
    param::Vector{QQMPolyRingElem}

    function CurveRationalParametrization(
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
    rat_param::Union{RationalParametrization, CurveRationalParametrization}

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

# A base point is a list of (LinearForm, Value) pairs
# defining the successive fiber cuts.
const BasePoint{T} = Vector{Tuple{QQMPolyRingElem, QQFieldElem}}

mutable struct RMnode{T <: QQMPolyRingElem}
    base_pt::BasePoint{T}
    polar_eqs::Vector{T}
    children::Vector{RMnode{T}}
end

mutable struct Roadmap{T <: QQMPolyRingElem}
    initial_ideal::Ideal{T}
    root::RMnode{T}
end

function _collect_roadmap(RMn::RMnode, F)
    data = [F(RMn)]
    for child in RMn.children
        append!(data, _collect_roadmap(child, F))
    end
    return data
end

function _fbr(I::Ideal{P}, base_pt::BasePoint{P}) where {P <: QQMPolyRingElem}
    @assert(!isempty(I.gens), "Empty polynomial vector")
    fb_eqs = [l - q for (l, q) in base_pt]
    return Ideal(vcat(I.gens, fb_eqs))
end

function all_eqs(RM::Roadmap)
    func(s) = _fbr(vcat(RM.initial_ideal.gens, s.polar_eqs) |> Ideal, s.base_pt)
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

mutable struct ParametricIdeal{K<:FracFieldElem}
    base_frac_field::FracField
    num_params::Int
    num_vars::Int
    symbols::Vector{Symbol}
    gens::Vector{MPoly{K}}
    sats::Vector{MPoly{K}}
    dim::Union{Int,Nothing}
    shift::Vector{Int}
    gb::Dict{UInt64,Vector{MPoly{K}}}
    gb_lock::ReentrantLock
    gb_spec::Dict{UInt64,Vector{QQMPolyRingElem}}
    gb_spec_lock::ReentrantLock
    mm::Dict{UInt64,MatSpaceElem{K}}
    mm_lock::ReentrantLock
    mm_spec::Dict{UInt64,QQMatrix}
    mm_spec_lock::ReentrantLock
    hm::Dict{UInt64,MatSpaceElem{K}}
    hm_lock::ReentrantLock
    hm_spec::Dict{UInt64,QQMatrix}
    hm_spec_lock::ReentrantLock
    prefer_gb::Bool
    radicalize::Bool
    gens_alt::Vector{MPoly{K}}

    # Generate magic numbers for shifting the interpolation points to avoid singularities
    # This works surprisingly well in practice, but it would be nice to have a more principled way to choose them
    function _magic_shift(n::Int)::Vector{Int}
        _magic_shift_cache = [11, 29, 47, 71, 97, 113, 149, 173, 197, 229]
        for i in length(_magic_shift_cache)+1:n
            push!(_magic_shift_cache, 2^i + rand(Int) % 2^i)
        end
        return [_magic_shift_cache[i] for i in 1:n]
    end

    function ParametricIdeal(gens::Vector{MPoly{K}};
        sats::Vector{MPoly{K}}=Vector{MPoly{K}}(),
        num_vars::Int=-1,
        prefer_gb::Bool=true,
        radicalize::Bool=false,
        gens_alt::Vector{MPoly{K}}=Vector{MPoly{K}}()
    ) where {K<:FracFieldElem}
        @assert !isempty(gens)
        I = new{K}()
        R = parent(gens[1])
        I.base_frac_field = base_ring(R)
        I.num_params = number_of_variables(base_ring(I.base_frac_field))
        I.num_vars = num_vars != -1 ? num_vars : number_of_variables(R)
        I.symbols = symbols(R)[end-I.num_vars+1:end]
        I.gens = gens
        I.sats = sats
        I.dim = nothing
        I.shift = _magic_shift(I.num_params)
        I.gb = Dict()
        I.gb_lock = ReentrantLock()
        I.gb_spec = Dict()
        I.gb_spec_lock = ReentrantLock()
        I.mm = Dict()
        I.mm_lock = ReentrantLock()
        I.mm_spec = Dict()
        I.mm_spec_lock = ReentrantLock()
        I.hm = Dict()
        I.hm_lock = ReentrantLock()
        I.hm_spec = Dict()
        I.hm_spec_lock = ReentrantLock()
        I.prefer_gb = prefer_gb
        I.radicalize = radicalize
        I.gens_alt = gens_alt
        return I
    end
end
