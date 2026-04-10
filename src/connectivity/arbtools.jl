"""
    arb_to_rat(x::ArbFieldElem) -> Tuple{Rational, Rational}

Return a rational interval enclosing the Arb ball.
"""
function arb_to_rat(x::ArbFieldElem)
    r = radius(x)
    lo = x - 2r
    hi = x + 2r
    return (simplest_rational_inside(lo), simplest_rational_inside(hi))
end


"""
    rat_to_arb(interval::Tuple, field::ArbField) -> ArbFieldElem

Convert rational interval (x1, x2) to Arb ball.
"""
function rat_to_arb(interval::Vector{QQFieldElem}, field::ArbField)
    x1, x2 = interval
    @assert x1 <= x2 "Invalid interval: x1 > x2"

    mid = field((x1 + x2) / 2)
    rad = field((x2 - x1) / 2)

    return ball(mid, rad)
end


"""
    evaluate_arb(f, x::ArbFieldElem)

Evaluate polynomial f at x using Arb arithmetic.
"""
function evaluate_arb(f, x::ArbFieldElem)
    is_zero(f) && return zero(parent(x))
    cf = coefficients_of_univariate(f)
    return evalpoly(x, cf)
end


"""
    evaluate_arb(f, g, x::ArbFieldElem)

Evaluate rational function f/g at x.
"""
function evaluate_arb(f, g, x::ArbFieldElem)
    @assert !is_zero(g) "Denominator must be non-zero"

    is_zero(f) && return zero(parent(x))

    cf = coefficients_of_univariate(f)
    cg = coefficients_of_univariate(g)

    num = evalpoly(x, cf)
    den = evalpoly(x, cg)

    return num / den
end


"""
    arb_eval(f, B::Vector{Tuple}, prec::Int)

Evaluate polynomial f over interval box B.
"""
function arb_eval(f, B::Vector{Tuple}, field::ArbField)
    BA = Vector{ArbFieldElem}(undef, length(B))
    for i in eachindex(B)
        BA[i] = rat_to_arb(B[i], field)
    end

    f_arb = change_base_ring(field, f)
    return evaluate(f_arb, BA)
end