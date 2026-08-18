import Nemo: identity_matrix, matrix

@testset "Algorithms -> Multiplication matrix" begin
    R, (x,) = polynomial_ring(QQ, [:x])
    f = [x^2 - 2 * x + 1]
    g = AlgebraicSolving.groebner_basis(AlgebraicSolving.Ideal(f))
    b = AlgebraicSolving._monomial_basis(g)
    M₁ = multiplication_matrix(g, b, one(R))
    @test M₁ == identity_matrix(QQ, 2)
    Mₓ = multiplication_matrix(g, b, x)
    @test Mₓ == matrix(QQ, [0 -1; 1 2])
end

@testset "Algorithms -> Parametric multiplication matrix" begin
    M = multiplication_matrix
    R, (a, b, c) = polynomial_ring(QQ, [:a, :b, :c])
    R′, (x,) = polynomial_ring(fraction_field(R), [:x], internal_ordering=:degrevlex)
    f = ParametricIdeal([a * x^2 + b * x + c])
    @test M(f, one(R′)) == matrix(fraction_field(R),
        [1 0; 0 1])
    @test M(f, x) == matrix(fraction_field(R),
        [0 -c//a; 1 -b//a])
    @test M(f, x^2) == matrix(fraction_field(R),
        [-c//a (b*c)//a^2; -b//a (b^2-a*c)//a^2])
end
