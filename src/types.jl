export Ideal

mutable struct RationalParametrization
    vars::Vector{Symbol}
    lf_cfs::Vector{fmpz}
    elim::fmpq_poly
    denom::fmpq_poly
    param::Vector{PolyElem}

    function RationalParametrization(
            vars::Vector{Symbol},
            lf_cfs::Vector{fmpz},
            elim::fmpq_poly,
            denom::fmpq_poly,
            param::Vector{PolyElem}
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

mutable struct Ideal{T <: MPolyElem}
    gens::Vector{T}
    dim::Int
    gb::Dict{Int, Vector{T}}
    real_sols::Vector{Vector{fmpq}}
    rat_param::RationalParametrization

    function Ideal(F::Vector{T}) where {T <: MPolyElem}
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
