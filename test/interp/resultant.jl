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

@testset "Interpolation -> Resultant and discriminant (trivariate)" begin
    I = AlgebraicSolving.Interpolation
    R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
    f = (x - y) * z^2 + (x * y + 1) * z + y^2 - 3
    g = (y + 1) * z^2 - x * z + 2 * x * y - x + 1
    for i in 1:3
        @test I.resultant(f, g, i) == I._flint_resultant(f, g, i)
        @test I.resultant(f, g, i; nr_thrds=4) == I._flint_resultant(f, g, i)
        @test I.discriminant(f, i) == I._flint_discriminant(f, i)
    end
end

@testset "Interpolation -> Resultant and discriminant (conventions, edge cases)" begin
    I = AlgebraicSolving.Interpolation
    R, (x, y) = polynomial_ring(QQ, [:x, :y])
    # degrees 1 to 4 in y (sign convention of the discriminant), nonconstant leading coefficients
    for f in [3y + x, x * y^2 + y - 7, -y^3 + QQ(1, 5) * x * y + 1, (x + 1) * y^4 - QQ(1, 3) * x * y^2 + 2]
        @test I.discriminant(f, 2) == I._flint_discriminant(f, 2)
        @test I.discriminant(f, 1) == I._flint_discriminant(f, 1)
    end
    # one polynomial constant in the eliminated variable
    @test I.resultant(x * y^2 + 1, x + 3, 2) == I._flint_resultant(x * y^2 + 1, x + 3, 2)
    @test I.resultant(x + 3, x * y^2 + 1, 2) == I._flint_resultant(x + 3, x * y^2 + 1, 2)
    # rational coefficients
    f = QQ(1, 3) * x^2 * y - QQ(1, 2) * y^2 + QQ(1, 7) * x
    g = QQ(5, 4) * x * y^2 + x^3 - 1
    @test I.resultant(f, g, 2) == I._flint_resultant(f, g, 2)
    @test I.resultant(f, g, 1) == I._flint_resultant(f, g, 1)
end

import Nemo: gens, coeff, constant_coefficient

@testset "Interpolation -> FLINT evaluation wrappers" begin
    I = AlgebraicSolving.Interpolation
    R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
    Qt, t = polynomial_ring(QQ, :t; cached=false)
    f = QQ(2, 3) * x^3 * y * z^2 - 5 * x * y^2 + QQ(1, 7) * z^4 - y + 3
    for (idx, vals) in [([1], [QQ(2)]), ([1, 3], [QQ(-1, 2), QQ(5)]), ([3, 2], [QQ(4), QQ(0)]), (Int[], QQFieldElem[])]
        @test I._flint_evaluate(f, idx, vals) == evaluate(f, gens(R)[idx], vals)
    end
    for i in 1:3
        others = [1:i-1; i+1:3]
        vals = [QQ(3, 4), QQ(-2)]
        g = evaluate(f, gens(R)[others], vals)
        @test I._flint_specialize(Qt, f, i, vals) == sum(constant_coefficient(coeff(g, [i], [k])) * t^k for k in 0:degree(g, i))
    end
    # mapping variables between rings by index (and to zero)
    S, (a, b) = polynomial_ring(QQ, [:a, :b])
    p = QQ(2, 3) * x^2 * z + 5 * x * y - z^3 + 1
    @test I._flint_map_gens(S, p, [1, 0, 2]) == QQ(2, 3) * a^2 * b - b^3 + 1
    @test I._flint_map_gens(R, QQ(2, 3) * a^2 * b - b^3 + 1, [1, 3]) == QQ(2, 3) * x^2 * z - z^3 + 1
    @test I._flint_map_gens(S, p, [2, 0, 1]) == QQ(2, 3) * a * b^2 - a^3 + 1
    S0, _ = polynomial_ring(QQ, Symbol[])
    @test I._flint_map_gens(S0, 2 * y + 7, [0, 0, 0]) == S0(7)
    # univariate input: the interpolation happens in a ring without variables
    R1, (u,) = polynomial_ring(QQ, [:u])
    @test I.resultant(u^3 - 2u + 1, 3u^2 + 5, 1) == I._flint_resultant(u^3 - 2u + 1, 3u^2 + 5, 1)
    @test I.discriminant(u^3 - 2u + 1, 1) == I._flint_discriminant(u^3 - 2u + 1, 1)
end
