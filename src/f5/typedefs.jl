# type for signature indices
const SigIndex = UInt16
# type for exponents in monomials
# ints because sigratios may require negative exponents
const Exp = UInt16
# types for hashvalue, ht index and divisor mask of a monomial
const MonIdx = Int32
const MonHash = UInt32
const DivMask = UInt32

struct Monomial{N}
    deg::Exp
    exps::SVector{N, Exp}
end
const Sig{N} = Tuple{SigIndex, Monomial{N}}

mutable struct LoadVector{P}
    elems::Vector{P}
    load::Int
end

abstract type Basis end

mutable struct POTBasis{N, C}<:Basis
    curr_indx::SigIndex
    sigs::Vector{Sig{N}}
    # elements in index k start in with index_cutoffs[k]
    index_cutoffs::Vector{Int}

    sigmasks::Vector{DivMask}

    lm_masks::Vector{DivMask}

    monomials::Vector{Vector{MonIdx}}
    coefficients::Vector{Vector{C}}

    syz_sigs::Vector{Monomial{N}}
    syz_masks::Vector{DivMask}
    # where do syzygies in current index start
    syz_curr_indx_start::Int
    syz_index_cutoffs::Vector{Int}

    basis_load::Int
    syz_load::Int
end

mutable struct SPair{N}
    top_sig::Monomial{N}
    bot_sig::Sig{N}
    top_sig_mask::DivMask
    bot_sig_mask::DivMask
    # top index = 0 -> pair is redundant
    top_index::Int
    bot_index::Int
end

const Pairset{N} = LoadVector{SPair{N}}

const ColIdx = UInt32
mutable struct MacaulayMatrix{C}

    # rows from upper, AB part of the matrix,
    # stored as vectors of corresponding exponents (already hashed and sorted)
    rows::Vector{Vector{ColIdx}}

    # maps column idx {1 ... ncols} to monomial hash idx {2 ... ht.load}
    # in some hashtable
    col2hash::Vector{MonIdx}

    # row coefficients
    coeffs::Vector{Vector{T}}

    #= sizes info =#
    # total number of allocated rows
    size::Int
    # number of pivots,
    # ie new basis elements discovered after matrix reduction
    npivots::Int
    # number of filled rows, nrows <= size
    nrows::Int
    # number of columns
    ncols::Int
end
