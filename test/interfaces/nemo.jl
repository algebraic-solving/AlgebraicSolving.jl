@testset "Interfaces -> Nemo" begin
    R, (x,y,z) = polynomial_ring(QQ,["x","y","z"], ordering=:degrevlex)
    F = [x^2+1-3, x*y-z, x*z^2-3*y^2]
    cmp = (Int32[2, 2, 2], BigInt[1, 1, -2, 1, 1, 1, -1, 1, 1, 1, -3, 1], Int32[2, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 2, 0, 2, 0])
    @test AlgebraicSolving._convert_to_msolve(F) == cmp
    R, (x,y,z) = polynomial_ring(GF(2147483659),["x","y","z"], ordering=:degrevlex)
    F = [x^2+1-3, x*y-z, x*z^2-3*y^2]
    # prime is bigger than 2^31, should throw an error
    @test_throws ErrorException AlgebraicSolving._convert_to_msolve(F)
    R, (x,y,z) = polynomial_ring(GF(101),["x","y","z"], ordering=:degrevlex)
    F = [x^2+1-3, x*y-z, x*z^2-3*y^2]
    res = AlgebraicSolving._convert_to_msolve(F)
    @test AlgebraicSolving._convert_finite_field_gb_to_abstract_algebra(Int32(3), res..., R) == F
end
