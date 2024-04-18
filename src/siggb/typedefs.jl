# type for signature indices
const SigIndex = UInt16
# type for exponents in monomials
# ints because sigratios may require negative exponents
const Exp = Int16
# types for hashvalue, ht index and divisor mask of a monomial
const MonIdx = Int32
const MonHash = UInt32
const DivMask = UInt32
# stuff for matrix
const ColIdx = UInt32
const Coeff = UInt32
# 64 bit buffer
const Cbuf = UInt64

struct Monomial{N}
    deg::Exp
    exps::SVector{N, Exp}
end
function monomial(exps::SV) where {N, SV <: StaticArray{Tuple{N}}}
    return Monomial{N}(sum(exps), exps)
end
const Sig{N} = Tuple{SigIndex, Monomial{N}}
const MaskSig = Tuple{SigIndex, DivMask}

const Polynomial = Tuple{Vector{Coeff}, Vector{MonIdx}}
const ModCache{N} = Dict{Tuple{Sig{N}, SigIndex}, Polynomial}

mutable struct Basis{N}
    sigs::Vector{Sig{N}}
    sigmasks::Vector{MaskSig}

    sigratios::Vector{Monomial{N}}

    # tree structure:
    #   - length of data for i is rewrite_nodes[i][1]
    #   - parent of i is rewrite_nodes[i][2]
    #   - children of i are rewrite_nodes[i][3:end]
    # careful: indices are shifted by + 1 compared to basis indices
    rewrite_nodes::Vector{Vector{Int}}

    lm_masks::Vector{DivMask}

    monomials::Vector{Vector{MonIdx}}
    coefficients::Vector{Vector{Coeff}}

    is_red::Vector{Bool}

    mod_rep_known::Vector{Vector{Bool}}
    mod_reps::Vector{Vector{Polynomial}}

    syz_sigs::Vector{Monomial{N}}
    syz_masks::Vector{MaskSig}

    # degrees of input elements
    degs::Vector{Exp}

    basis_load::Int
    basis_size::Int
    # for storing the initial polynomials
    input_load::Int
    input_size::Int
    basis_offset::Int
    syz_load::Int
    syz_size::Int
end

# TODO: should these be stored in a more vectorized way?
mutable struct SPair{N}
    top_sig::Sig{N}
    bot_sig::Sig{N}
    top_sig_mask::DivMask
    bot_sig_mask::DivMask
    # top index = 0 -> pair is redundant
    top_index::Int
    bot_index::Int
    deg::Exp
end

mutable struct Pairset{N}
    elems::Vector{SPair{N}}
    load::Int
    size::Int
end

mutable struct MacaulayMatrix{N}

    # stored as vectors of corresponding exponents (already hashed and sorted)
    rows::Vector{Vector{MonIdx}}

    # pivot row index for colidx is pivots[colidx] 
    pivots::Vector{Int}
    pivot_size::Int

    sigs::Vector{Sig{N}}
    parent_inds::Vector{Int}
    # sig(row[i]) < sig(row[j]) <=> sig_order[i] < sig_order[j]
    sig_order::Vector{Int}

    # maps column index to corresponding hash index
    col2hash::Vector{ColIdx}

    # row coefficients
    coeffs::Vector{Vector{Coeff}}

    #= sizes info =#
    # total number of allocated rows
    size::Int
    # number of filled rows, nrows <= size
    nrows::Int
    # number of columns
    ncols::Int

    # for each i in toadd rows[i] should be added to basis/syzygies
    toadd::Vector{Int}
    toadd_length::Int
end

# struct to remember the row reductions we did
mutable struct TracerMatrix
    # first index row index, second one rewr ind
    rows::Dict{Sig, Tuple{Int, Int}}
    # if row i has been added to basis
    is_basis_row::Dict{Int, Int}
    row_ind_to_sig::Dict{Int, Sig}
    diagonal::Vector{Coeff}
    col_inds_and_coeffs::Vector{Vector{Tuple{Int, Coeff}}}
end

mutable struct Tracer
    mats::Vector{TracerMatrix}
    basis_ind_to_mat::Vector{Int}
    syz_ind_to_mat::Vector{Int}
    load::Int
    size::Int
end

# For Index ordering
mutable struct IndOrder
    ord::Vector{SigIndex}
    incompat::Dict{Tuple{SigIndex, SigIndex}, Bool}
    max_ind::SigIndex
end

# This is to store where certain elements come from
const Tags = Dict{SigIndex, Symbol}

# For syzygy processing during splitting
const SyzInfo = Tuple{SigIndex, Dict{SigIndex, Bool}}

# For output of decomp algorithms
abstract type AffineCell end

mutable struct LocClosedSet{T<:MPolyRingElem} <: AffineCell
    seq::Vector{T}
    ineqns::Vector{Vector{T}}
    gbs::Vector{Vector{T}}

    function LocClosedSet{T}(seq::Vector{T}) where {T<:MPolyRingElem}
        @assert !isempty(seq) "cannot construct affine cell from no equations."
        R = parent(first(seq))
        gb = saturate(seq, last(gens(R)))
        return new(seq, [[last(gens(R))]], [gb])
    end
end

mutable struct WLocClosedSet{T<:MPolyRingElem} <: AffineCell
    seq::Vector{T}
    tags::Vector{Symbol}
    hull_eqns::Vector{T}
    hypplanes::Vector{T}
    codim_upper_bound::Int
    ineqns::Vector{T}
    gbs::Vector{Vector{T}}

    function WLocClosedSet{T}(seq::Vector{T}) where {T<:MPolyRingElem}
        @assert !isempty(seq) "cannot construct affine cell from no equations."
        R = parent(first(seq))
        codim_upper_bound = min(length(seq), ngens(R) - 1)
        hypplanes = [random_lin_comb(gens(R)) for _ in 1:(ngens(R) - codim_upper_bound)]
        gb = saturate(vcat(seq, hypplanes), last(gens(R)))
        return new(seq, [:seq for _ in 1:length(seq)],
                   T[], hypplanes, codim_upper_bound,
                   [last(gens(R))], [gb])
    end
end

# for benchmarking
mutable struct Timings
    sym_pp_time::Float32
    lin_alg_time::Float32
    select_time::Float32
    update_time::Float32
    module_time::Float32
    comp_lc_time::Float32
end
