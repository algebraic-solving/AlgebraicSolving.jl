@testset "Algorithms -> Gröbner bases" begin
    R, (x,y,z) = PolynomialRing(QQ,["x","y","z"], ordering=:degrevlex)
    F = [x^2+1-3, x*y-z, x*z^2-3*y^2]
    #= not a finite field =#
    @test_throws ErrorException groebner_basis(F)
    R, (x,y,z) = PolynomialRing(GF(101),["x","y","z"], ordering=:degrevlex)
    F = [x+2*y+2*z-1, x^2+2*y^2+2*z^2-x, 2*x*y+2*y*z-y]
    G = groebner_basis(F)
    H = MPolyElem[
         x + 2*y + 2*z + 100
         y*z + 82*z^2 + 10*y + 40*z
         y^2 + 60*z^2 + 20*y + 81*z
         z^3 + 28*z^2 + 64*y + 13*z
                 ]
    @test G == H
    G = groebner_basis(F, eliminate=2)
    H = MPolyElem[
         z^4 + 38*z^3 + 95*z^2 + 95*z
                 ]
    @test G == H
end
