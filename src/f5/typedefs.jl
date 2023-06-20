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

struct Monomial{N}
    deg::Exp
    exps::SVector{N, Exp}
end
function monomial(exps::SVector{N, Exp}) where N
    return Monomial{N, Exp}(sum(exps), exps)
end
const Sig{N} = Tuple{SigIndex, Monomial{N}}
const MaskSig = Tuple{SigIndex, DivMask}

mutable struct LoadVector{P}
    elems::Vector{P}
    load::Int
end

mutable struct Basis{N, C}
    sigs::Vector{Monomial{N}}
    sigmasks::Vector{MaskSig}

    sigratios::Vector{Monomial{N}}

    lm_masks::Vector{DivMask}

    monomials::Vector{Vector{MonIdx}}
    coefficients::Vector{Vector{C}}

    is_red::Vector{Bool}

    syz_sigs::Vector{Monomial{N}}
    syz_masks::Vector{MaskSig}

    # degrees of input elements
    # TODO: do we ever need this
    degs::Vector{Exp}

    basis_load::Int
    syz_load::Int
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

const Pairset{N} = LoadVector{SPair{N}}

mutable struct MacaulayMatrix{C}

    # stored as vectors of corresponding exponents (already hashed and sorted)
    rows::Vector{Vector{ColIdx}}

    # pivot row index for colidx is pivots[colidx] 
    pivots::Vector{Int}
    pivot_size::Int

    # sig(row[i]) < sig(row[j]) <=> sig_order[i] < sig_order[j]
    sig_order::Vector{Int}

    # maps column idx {1 ... ncols} to monomial hash idx {2 ... ht.load}
    # in some hashtable
    col2hash::Vector{MonIdx}

    # row coefficients
    coeffs::Vector{Vector{T}}

    #= sizes info =#
    # total number of allocated rows
    size::Int
    # number of filled rows, nrows <= size
    nrows::Int
end
