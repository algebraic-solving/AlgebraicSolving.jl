@testset "Interfaces -> Nemo" begin
    R, (x,y,z) = polynomial_ring(QQ,["x","y","z"], internal_ordering=:degrevlex)
    F = [x^2+1-3, x*y-z, x*z^2-3*y^2]
    cmp = (Int32[2, 2, 2], BigInt[1, 1, -2, 1, 1, 1, -1, 1, 1, 1, -3, 1], Int32[2, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 2, 0, 2, 0], 3)
    @test AlgebraicSolving._convert_to_msolve(F) == cmp
    # issue #54
    R, (x1, x2) = polynomial_ring(GF(17), ["x1", "x2"])
    I = Ideal([x1, R(0)])
    cmp = (Int32[1], BigInt[1], Int32[1, 0], 1) 
    @test AlgebraicSolving._convert_to_msolve(I.gens) == cmp
    for _GF in [GF, AlgebraicSolving.Nemo.Native.GF]
        R, (x,y,z) = polynomial_ring(GF(2147483659),["x","y","z"], internal_ordering=:degrevlex)
        F = [x^2+1-3, x*y-z, x*z^2-3*y^2]
        # prime is bigger than 2^31, should throw an error
        @test_throws ErrorException AlgebraicSolving._convert_to_msolve(F)
        R, (x,y,z) = polynomial_ring(GF(101),["x","y","z"], internal_ordering=:degrevlex)
        F = [x^2+1-3, x*y-z, x*z^2-3*y^2]
        res = AlgebraicSolving._convert_to_msolve(F)
        @test AlgebraicSolving._convert_finite_field_array_to_abstract_algebra(Int32(3), res[1], res[2], res[3], R) == F
    end
    R, (x,y,z) = polynomial_ring(QQ,["x","y","z"], internal_ordering=:degrevlex)
    F = [x^2+1-3, x*y-z, x*z^2-3*y^2]
    res = AlgebraicSolving._convert_to_msolve(F)
    res_qq =AlgebraicSolving.QQFieldElem[res[2][i] for i in 1:2:length(res[2])]
    @test AlgebraicSolving._convert_rational_array_to_abstract_algebra(Int32(3), res[1], res_qq, res[3], R) == F
end
