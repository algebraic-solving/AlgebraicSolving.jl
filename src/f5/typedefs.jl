# type for signature indices
const SigIndex = UInt16
# type for exponents in monomials
# ints because sigratios may require negative exponents
const Exp = Int16
# types for hashvalue and divisor mask of a monomial
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

mutable struct POTBasis{N, C}
    curr_indx::SigIndex
    sigs::Vector{Sig{N}}
    # elements in index k start in with index_cutoffs[k]
    index_cutoffs::Vector{Int}

    sigmasks::Vector{DivMask}

    lm_masks::Vector{DivMask}

    monomials::Vector{Vector{MonHash}}
    coefficients::Vector{Vector{C}}

    syz_sigs::Vector{Monomial{N}}
    syz_masks::Vector{DivMask}
    # where do syzygies in current index start
    syz_curr_indx_start::Int
    syz_index_cutoffs::Vector{Int}

    basis_load::Int
    syz_load::Int
end

mutable struct POTSPair{N}
    top_sig::Monomial{N}
    bot_sig::Sig{N}
    top_sig_mask::DivMask
    bot_sig_mask::DivMask
    # top index = 0 -> pair is redundant
    top_index::Int
    bot_index::Int
end

const POTPairset{N} = LoadVector{POTSPair{N}}

