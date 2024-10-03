# compute normal form with respect to basis
# *without* computing a GB for the corresponding ideal
function normal_form(F::Vector{T}, gb::Vector{T}) where {T <: MPolyRingElem}
    I = Ideal(gb)
    I.gb[0] = gb
    return normal_form(F, I)
end

function normal_form(mons::Vector{MonIdx},
                        coeffs::Vector{Coeff},
                        ht::MonomialHashtable{N},
                        gb::Vector{<:MPolyRingElem}) where N
    

    R = parent(first(gb))
    @inbounds mons = [ht.exponents[midx] for midx in mons]
    f = convert_to_pol(R, mons, coeffs)
    F = [f]
    return normal_form(F, gb)
end
