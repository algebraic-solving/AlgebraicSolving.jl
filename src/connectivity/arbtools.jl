"""
    arb_to_rat(x::ArbFieldElem) -> Tuple{Rational, Rational}

Return a rational interval enclosing the Arb ball.
"""
function arb_to_rat(x::ArbFieldElem)
    r = radius(x)
    for i in 2:5
        lo = x - i*r
        hi = x + i*r
        I = map(simplest_rational_inside, [lo, hi])
        contains(rat_to_arb(I, parent(x)), x) && return I
    end
    error("Problem in boxes computations with Arb. Try increasing precision.")
end


"""
    rat_to_arb(interval::Vector{QQFieldElem}, field::ArbField) -> ArbFieldElem

Convert rational interval (x1, x2) to Arb ball.
"""
function rat_to_arb(interval::Vector{QQFieldElem}, field::ArbField)
    x1, x2 = interval
    @assert x1 <= x2 "Invalid interval: x1 > x2"
    new_field = field
    for i in 2:5
        mid = field((x1 + x2) / 2)
        rad = field((x2 - x1) / 2)
        I = ball(mid, rad)
        all(contains(I, x) for x in interval) && return I
        new_field = ArbField(i * field.prec)
    end
    error("Problem in boxes computations with Arb. Try increasing precision.")
end


"""
    evaluate_arb(f, g, x::Vector{QQFieldElem})

Evaluate rational function f/g at x.
Return a rational interval enclosing the result, computed using Arb.
"""
function evaluate_arb(f, g, x::Vector{QQFieldElem}, field::ArbField)
    @assert !is_zero(g) "Denominator must be non-zero"

    is_zero(f) && return zero(parent(x))

    cf = field.(coefficients_of_univariate(f))
    cg = field.(coefficients_of_univariate(g))

    x_arb = rat_to_arb(x, field)
    num = evalpoly(x_arb, cf)
    den = evalpoly(x_arb, cg)

    return arb_to_rat(num / den)
end


"""
    arb_eval(f, B::Vector{Tuple}, prec::Int)

Evaluate polynomial f over interval box B.
Return an Arb ball containing the image of B by f.
"""
function arb_eval(f, B::Vector{Tuple}, field::ArbField)
    BA = Vector{ArbFieldElem}(undef, length(B))
    for i in eachindex(B)
        BA[i] = rat_to_arb(B[i], field)
    end

    f_arb = change_base_ring(field, f)
    return evaluate(f_arb, BA)
end