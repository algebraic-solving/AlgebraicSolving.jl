# functions for polynomials
function add_pols(mons1::Vector{M},
                  mons2::Vector{M},
                  coeffs1::Vector{Coeff},
                  coeffs2::Vector{Coeff},
                  vch::Val{Char}) where {M <: Monomial, Char}

    l1 = length(mons1)
    l2 = length(mons2)
    mons_res = Vector{M}(undef, l1 + l2)
    coeffs_res = Vector{Coeff}(undef, l1 + l2)
    
    ind1 = 1
    ind2 = 1
    new_l = 0 
    @inbounds while ind1 <= l1 && ind2 <= l2
        new_l += 1
        m1 = mons1[ind1]
        m2 = mons2[ind2]
        if m1 == m2
            mons_res[new_l] = m1
            coeffs_res[new_l] = add(coeffs1[ind1], coeffs2[ind2], vch)
            ind1 += 1
            ind2 += 1
        elseif lt_drl(m2, m1)
            mons_res[new_l] = m1
            coeffs_res[new_l] = coeffs1[ind1]
            ind1 += 1
        else
            mons_res[new_l] = m2
            coeffs_res[new_l] = coeffs2[ind2]
            ind2 += 1
        end
    end

    resize!(mons_res, new_l)
    resize!(coeffs_res, new_l)

    return mons_res, coeffs_res
end
