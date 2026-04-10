function isolate(f; prec = 32)
    # Univariate isolation with msolve
    total_degree(f) == 0 && return []
    @assert is_univariate(f) "Not univariate polynomial"

    sols = real_solutions(Ideal([change_ringvar(f)]), precision=prec, interval=true, info_level=0)
    return getindex.(sols, 1)
end

function isolate_eval(f, ivar, val; prec=64)
    # univariate isolation of roots of a bivariate polynomial f
    # whose ivar-th variable is evaluated at val
    fev = change_ringvar(evaluate(f, [ivar], [val]))
    #fev *= fev |> coefficients .|> denominator |> lcm
    return isolate(fev, prec=prec)
end