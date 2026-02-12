@testset "Interpolation -> Newton" begin
    R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
    n = 5
    xs = [ZZ(i) for i in 1:n]
    for f in [zero(R), one(R), x, y, z, x + z, y * z, x * y * (x + y + z)]
        @test newton(R, t -> evaluate(f, t), [n - 1, n - 1, n - 1]) == f
    end
end
