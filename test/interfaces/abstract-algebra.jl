@testset "Interfaces -> AbstractAlgebra" begin
    R, (x,y,z) = PolynomialRing(QQ,["x","y","z"], ordering=:degrevlex)
    F = [x^2+1-3, x*y-z, x*z^2-3*y^2]
    cmp = (Int32[2, 2, 2], BigInt[1, 1, -2, 1, 1, 1, -1, 1, 1, 1, -3, 1], Int32[2, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 2, 0, 2, 0])
    @test AlgebraicSolving.convert_to_msolve(F) == cmp
end
