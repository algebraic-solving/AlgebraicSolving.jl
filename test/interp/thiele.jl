@testset "Interpolation -> Thiele" begin
    R, x = polynomial_ring(QQ, :x)
    for f in [zero(R), one(R), x, x * x, x // (x + 1), (x^2 + 1) // (x^4 + x^3 + 1)]
        @test thiele(R, t -> evaluate(f, t)) == f
    end
end
