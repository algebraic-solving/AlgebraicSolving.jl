@testset "Interpolation -> Resultant" begin
    R, (x, y) = polynomial_ring(QQ, [:x, :y])

    f = 3x * y^2 + (x + 1) * y + 3
    g = 6(x + 1) * y + (x^3 + 2x + 2)

    @test AlgebraicSolving.Interpolation.resultant(f, g, 2) == 3 * x^7 + 6 * x^5 - 6 * x^3 + 96 * x^2 + 192 * x + 96
end

@testset "Interpolation -> Discriminant" begin
    R, (x, y) = polynomial_ring(QQ, [:x, :y])

    f = x * y^2 + (x + 1) * y + 3

    @test AlgebraicSolving.Interpolation.discriminant(f, 2) == x^2 - 10 * x + 1
end
