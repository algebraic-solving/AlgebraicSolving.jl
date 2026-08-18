import Nemo: fraction_field

@testset "Algorithms -> Parametric Gröbner basis" begin
    R, (a, b, c) = polynomial_ring(QQ, [:a, :b, :c])
    R′, (x,) = polynomial_ring(fraction_field(R), [:x], internal_ordering=:degrevlex)
    f = ParametricIdeal([a * x^2 + b * x + c])
    gb = groebner_basis(f)
    @test gb == [(a * x^2 + b * x + c) / a]
end

@testset "Algorithms -> Monomial basis" begin
    R, (x,) = polynomial_ring(QQ, [:x])
    f = [x^2 - 2 * x + 1]
    g = AlgebraicSolving.groebner_basis(AlgebraicSolving.Ideal(f))
    b = AlgebraicSolving._monomial_basis(g)
    @test length(b) == 2
end
