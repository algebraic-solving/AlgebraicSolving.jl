include("../siggb/helpers.jl")
include("dimension.jl")

export affine_hilbert_series, hilbert_series

function _hilbert_series_mono(exps::Vector{Vector{Int}})
    r = sort!(exps, by=reverse) |> length
    println(exps)
    A, t = polynomial_ring(ZZ, 't')
    h = iszero(r) ? 1 : 1-t^(sum(exps[1]))
    for i in 2:r
        # Compute generators for (x^a1,...,x^a{i-1}):x^ai
        sat = [ [ max(a[j]-exps[i][j], 0) for j in eachindex(a)]
                                                for a in exps[1:i-1] ]
        # Reduce to minimal generators
        sat = [sat[j] for j in eachindex(sat) if 
                      !any(all(sat[k] .<= sat[j]) for k in eachindex(sat) if k!=j)]

        # TODO: why not removing pure powers xi^k?
        nolin_sat = [u for u in sat if sum(u)>1]
        hsat = (1-t)^(length(sat)-length(nolin_sat))*_hilbert_series_mono(nolin_sat)
        h = 1-t^(sum(exps[i]))*hsat
    end
    return h
end

function hilbert_series(I)
    gb = get(I.gb, 0, groebner_basis(I, complete_reduction = true))
    lexps = (_drl_lead_exp).(gb)
    return _hilbert_series_mono(lexps)
end

function affine_hilbert_series(I)
    # TODO: homogenize gb if exists
    F = I.gens
    Fh = homogenize(F)
    return hilbert_series(Ideal(Fh))
end