@testset "Algorithms -> Gröbner bases" begin
    R, (x,y,z) = polynomial_ring(QQ,["x","y","z"], ordering=:degrevlex)
    F = [x^2+1-3, x*y-z, x*z^2-3*y^2]
    #= not a finite field =#
    @test_throws ErrorException groebner_basis(Ideal(F))
    R, (x,y,z) = polynomial_ring(GF(101),["x","y","z"], ordering=:degrevlex)
    I = Ideal([x+2*y+2*z-1, x^2+2*y^2+2*z^2-x, 2*x*y+2*y*z-y])
    G = groebner_basis(I)
    H = MPolyRingElem[
         x + 2*y + 2*z + 100
         y*z + 82*z^2 + 10*y + 40*z
         y^2 + 60*z^2 + 20*y + 81*z
         z^3 + 28*z^2 + 64*y + 13*z
                 ]
    @test G == H
    @test I.gb[0] == H

    G = groebner_basis(I, eliminate=2)
    H = MPolyRingElem[
         z^4 + 38*z^3 + 95*z^2 + 95*z
        ]
    @test G == H

    @test_throws ErrorException eliminate(I,0)
    L = eliminate(I,2)
    @test L == H
    @test I.gb[2] == H
end

@testset "Algorithms -> Sig Gröbner bases" begin
    R, (x,y,z) = polynomial_ring(QQ,["x","y","z"], ordering=:degrevlex)
    F = [x^2+1-3, x*y-z, x*z^2-3*y^2]
    #= not a finite field =#
    @test_throws ErrorException sig_groebner_basis(F)
    R, (x,y,z) = polynomial_ring(GF(17),["x","y","z"], ordering=:degrevlex)
    F = [x^2+1-3, x*y-z, x*z^2-3*y^2]
    #= not homogeneous =#
    @test_throws ErrorException sig_groebner_basis(F)

    #= GB test 1 =#
    Fhom = AlgebraicSolving._homogenize(F)
    sgb = sig_groebner_basis(Fhom)
    @test AlgebraicSolving._is_gb(sgb)

    #= GB test 2 =#
    R, (x,y,z,w) = polynomial_ring(GF(65521),["x","y","z","w"], ordering=:degrevlex)
    F = cyclic(R).gens
    Fhom = AlgebraicSolving._homogenize(F)
    sgb = sig_groebner_basis(Fhom)
    @test AlgebraicSolving._is_gb(sgb)

    #= GB test 3 =#
    F = katsura(R).gens
    Fhom = AlgebraicSolving._homogenize(F)
    sgb = sig_groebner_basis(Fhom)
    @test AlgebraicSolving._is_gb(sgb)
end
