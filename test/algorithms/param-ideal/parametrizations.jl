import Nemo.Generic: FracFieldElem, MPoly

@testset "Algorithms -> Parametric rational parametrization" begin
    R, (a, b, c) = polynomial_ring(QQ, [:a, :b, :c])
    R′, (x,) = polynomial_ring(fraction_field(R), [:x], internal_ordering=:degrevlex)
    f = ParametricIdeal([a * x^2 + b * x + c])
    @test rational_parametrization(f) == (
        [:x],
        ZZRingElem[],
        x^2 + b // a * x + c // a,
        x + (1 // 2 * b) // a,
        MPoly{FracFieldElem{QQMPolyRingElem}}[]
    )
end

@testset "Algorithms -> Parametric rational parametrization (unlucky points)" begin
    R, (a, b, c) = polynomial_ring(QQ, [:a, :b, :c])
    R′, (x,) = polynomial_ring(fraction_field(R), [:x], internal_ordering=:degrevlex)
    for shift in ([-1, 2, 2], [-2, 1, 1])
        f = ParametricIdeal([a * x^2 + b * x + c])
        f.shift = shift
        @test rational_parametrization(f) == (
            [:x],
            ZZRingElem[],
            x^2 + b // a * x + c // a,
            x + (1 // 2 * b) // a,
            MPoly{FracFieldElem{QQMPolyRingElem}}[]
        )
    end
end
