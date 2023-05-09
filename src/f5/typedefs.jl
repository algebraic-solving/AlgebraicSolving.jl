const SigIndex = UInt16
const MonIndex = Int32
# store degree in first exponent
const MonHash = UInt32
const Sig{N} = Tuple{SigIndex, Monomial{N}}

include("monomials_new.jl")

# for readibility
index(a::Sig) = a[1]
monomial(a::Sig) = a[2]

# TODO: should we assume that things are sorted by sig degree
mutable struct POTBasis{N, C}
    curr_indx::SigIndex
    # where basis elements in current index start
    sigs::Vector{Sig{N}}
    curr_indx_start::Int
    # elements in index i start in index_cutoffs[i][1] and end in index_cutoffs[i][2]
    index_cutoffs::Dict{SigIndex, Tuple{Int, Int}}

    sigmasks::Vector{DivMask}
    sigratios::Vector{Monomial{N}}

    # TODO: should we store lms outside the hashtable?
    lms::Vector{MonHash}
    lm_masks::Vector{DivMask}

    monomials::Vector{Vector{MonHash}}
    coefficients::Vector{Vector{C}}

    syz_sigs::Vector{Monomial{N}}
    syz_masks::Vector{DivMask}
    # where do syzygies in current index start
    syz_curr_indx_start::Int
    syz_index_cutoffs::Dict{SigIndex, Tuple{Int, Int}}

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

mutable struct Pairset{P}
    pairs::Vector{P}
    load::Int
end

function update_pairset!(pairset::Pairset{POTSPair{N}},
                         basis::POTBasis,
                         basis_ht::MonomialHashtable,
                         new_basis_idx::Int) where N


    # resize pairset if needed
    pair_size = length(pairset.pairs)
    num_new_pairs = new_basis_idx - 1
    if pairset.load + num_new_pairs >= pair_size
        resize!(pairset.pairs, max(2 * pair_size, pair.load - num_new_pairs))
    end

    new_sig_mon = monomial(basis.sigs[new_basis_idx])

    # check existing pairs for rewriteability
    bmask = divmask(new_sig_mon, basis_ht.divmap, basis_ht.ndivbits)
    @inbounds for i in 1:(pairset.load)
        p = pairset.pairs[i]
        iszero(p.top_index) && continue
        if mask_div(bmask, p.top_sig_mask)
            if div(new_sig_mon, p.top_sig)
                pairset.pairs[i].top_index = 0
                continue
            end
        end
        basis.curr_indx != index(p.bot_sig) && continue
        if mask_div(bmask, p.bot_sig_mask)
            if div(new_sig_mon, monomial(p.bot_sig))
                pairset.pairs[i].top_index = 0
                continue
            end
        end
    end

    new_lm = basis_ht.exponents[basis.lms[new_basis_idx]]
    # pair construction loop
    # TODO: where to check that S-pair is regular?
    @inbounds for i in 1:(new_basis_idx - 1)
        basis_lm = basis_ht.exponents[basis.lms[i]]
        # lcm(new_lm, basis_lm)/new_lm
        mult_new_elem = lcm_div(new_lm, basis_lm)
        new_pair_sig = mul(mult_new_elem, new_sig_mon)

        # check multiplied signature of new element against syzygies
        # TODO: join the syzygy checks into one loop
        rewr_syz(new_pair_sig, basis, basis.syz_curr_indx_start, basis.syz_load,
                 basis.curr_indx_start - 1) && continue

        mult_basis_elem = lcm_div(basis_lm, new_lm)
        basis_pair_sig_mon = mul(mult_basis_elem, monomial(basis.sigs[i]))
        ind = index(basis.sigs[i])

        ind_cuts = basis.index_cutoffs[ind]
        # check multiplied signature of basis element against syzygies
        rewr_syz(basis_pair_sig_mon, basis.syz_index_cutoffs[ind]...,
                 ind_cuts[1] - 1) && continue

        basis_pair_mon_mask = divmask(basis_pair_sig_mon, basis_ht.divmap,
                                     basis_ht.ndivbits)
        is_rewr = false
        # check multiplied signature of basis element against basis sigs
        @inbounds for j in ind_cuts[1]:i-1
            if mask_div(basis.sigmasks[j], basis_pair_mon_mask)
                is_rewr = div(monomial(basis.sigs[j]), basis_pair_sig_mon)
                is_rewr && break
            end
        end
        is_rewr && continue
        
        new_sigratio = basis.sigratios[new_basis_idx]
        basis_sigratio = basis.sigratios[i]
        new_pair_sig_mask = divmask(new_pair_sig, basis_ht.divmap,
                                    basis_ht.ndivbits)
        new_pair = if index(basis.sigs[i]) < basis.curr_indx || lt_drl(basis_sigratio, new_sigratio)
            POTSPair(new_pair_sig, Sig(index(basis.sigs[i]), basis_pair_sig_mon),
                     new_pair_sig_mask, basis_pair_mon_mask, new_basis_idx, i)
        else
            POTSPair(basis_pair_sig_mon, new_pair_sig,
                     basis_pair_mon_mask, new_pair_sig_mask, i, new_basis_idx)
        end
            
        pairset.pairs[pairset.load + 1] = new_pair
        pairset.load += 1
    end
end

function rewr_syz(sig_mon::Monomial{N},
                  basis::POTBasis{N},
                  syz_indx_start,
                  syz_indx_end,
                  basis_indx_end) where N

    sig_mon_mask = divmask(sig1_mon, basis_ht.divmap,
                                basis_ht.ndivbits)

    # rewrite check against syzygies
    @inbounds for j in syz_indx_start:syz_indx_end
        if mask_div(basis.syz_masks[j], sig_mon_mask)
            div(monomial(basis.syz_sigs[j]), sig_mon) && return true
        end
    end

    # rewrite check against koszul syzygies
    @inbounds for j in 1:basis_indx_end
        if mask_div(basis.lm_masks[j], sig_mon_mask)
            div(basis_ht.exponents[basis.lms[j]], sig_mon) && return true
        end
    end

   return false 
end
