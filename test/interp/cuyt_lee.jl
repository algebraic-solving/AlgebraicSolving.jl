import Nemo: coeff

@testset "Interpolation -> Cuyt-Lee" begin
    R, (z1, z2) = polynomial_ring(QQ, [:z1, :z2])
    f = (3 + 2 * z1 + 4 * z2 + 7 * z1^2 + 5 * z1 * z2 + 6 * z2^2) // (1 + 7 * z1 + 8 * z2 + 10 * z1^2 + z1 * z2 + 9 * z2^2)
    @test cuyt_lee(R, t -> evaluate(f, t)) == f
    f = (3 + 2 * z1 + 4 * z2 + 7 * z1^2 + 5 * z1 * z2 + 6 * z2^2) // (z1 + z2 + 10 * z1^2 + z1 * z2 + 9 * z2^2)
    @test cuyt_lee(R, t -> evaluate(f, t)) == f
end

@testset "Interpolation -> Normal forms via Cuyt-Lee" begin
    R, (a, b, c) = polynomial_ring(QQ, [:a, :b, :c])
    _, (x,) = polynomial_ring(QQ, [:x])
    @test cuyt_lee(R, t -> coeff(divrem(x^3, [t[1] * x^2 + t[2] * x + t[3]])[2], 1)) == (b^2 - a * c) // a^2
end
