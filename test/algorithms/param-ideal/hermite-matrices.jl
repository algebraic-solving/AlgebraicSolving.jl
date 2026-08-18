import Nemo: zero_matrix

@testset "Algorithms -> Hermite matrix 1" begin
    R, (x,) = polynomial_ring(QQ, [:x])
    f = [R(1)]
    g = AlgebraicSolving.groebner_basis(AlgebraicSolving.Ideal(f))
    b = AlgebraicSolving._monomial_basis(g)
    H₁ = hermite_matrix(g, b, one(R))
    @test H₁ == zero_matrix(QQ, 0, 0)
end

@testset "Algorithms -> Hermite matrix 2" begin
    R, (x,) = polynomial_ring(QQ, [:x])
    f = [x^2 - 2 * x + 1]
    g = AlgebraicSolving.groebner_basis(AlgebraicSolving.Ideal(f))
    b = AlgebraicSolving._monomial_basis(g)
    H₁ = hermite_matrix(g, b, one(R))
    @test H₁ == matrix(QQ, [2 2; 2 2])
    Hₓ = hermite_matrix(g, b, x)
    Mₓ = multiplication_matrix(g, b, x)
    @test Hₓ == H₁ * Mₓ
end

@testset "Algorithms -> Parametric Hermite matrix" begin
    R, (a, b, c) = polynomial_ring(QQ, [:a, :b, :c])
    R′, (x,) = polynomial_ring(fraction_field(R), [:x], internal_ordering=:degrevlex)
    f = ParametricIdeal([a * x^2 + b * x + c])
    g = one(R′)
    H = hermite_matrix(f, g)
    @test H == matrix(fraction_field(R), [2 -b//a; -b//a (b^2-2*a*c)//a^2])
end

@testset "Algorithms -> Parametric Hermite matrix with zero parameters" begin
    R, () = polynomial_ring(QQ, Symbol[])
    R′, (x,) = polynomial_ring(fraction_field(R), [:x], internal_ordering=:degrevlex)
    f = ParametricIdeal([x^2 - 2 * x + 1])
    g = one(R′)
    H₁ = hermite_matrix(f, g)
    @test H₁ == matrix(QQ, [2 2; 2 2])
end
