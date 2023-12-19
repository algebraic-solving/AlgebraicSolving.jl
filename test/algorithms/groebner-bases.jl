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

    I = Ideal([R(0)])
    G = groebner_basis(I)
    @test G == [R(0)]
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

    #= GB test 4 (pivot setting bug) =#
    R, (x1, x2, x3, x4) = polynomial_ring(GF(65521), ["x1", "x2", "x3", "x4"], ordering=:degrevlex)
    F = [11523*x1^4 + 30378*x1^3*x2 + 30154*x1^2*x2^2 + 10157*x1*x2^3 - 28136*x2^4 - 4771*x1^3*x3 - 21056*x1^2*x2*x3 + 15696*x1*x2^2*x3 - 16144*x2^3*x3 - 1553*x1^2*x3^2 - 30379*x1*x2*x3^2 - 12735*x2^2*x3^2 + 18058*x1*x3^3 + 24670*x2*x3^3 - 16379*x3^4 + 24196*x1^3*x4 - 19411*x1^2*x2*x4 + 17610*x1*x2^2*x4 - 5715*x2^3*x4 - 21186*x1^2*x3*x4 - 22865*x1*x2*x3*x4 - 1939*x2^2*x3*x4 - 5685*x1*x3^2*x4 + 8508*x2*x3^2*x4 + 21819*x3^3*x4 - 24868*x1^2*x4^2 - 18233*x1*x2*x4^2 - 14116*x2^2*x4^2 + 28291*x1*x3*x4^2 - 9068*x2*x3*x4^2 - 15138*x3^2*x4^2 + 8921*x1*x4^3 - 18808*x2*x4^3 - 3005*x3*x4^3 + 7368*x4^4,
         31703*x1^4 + 23616*x1^3*x2 + 20696*x1^2*x2^2 - 7125*x1*x2^3 + 15334*x2^4 + 26619*x1^3*x3 + 2173*x1^2*x2*x3 - 31312*x1*x2^2*x3 - 31386*x2^3*x3 - 25244*x1^2*x3^2 - 28729*x1*x2*x3^2 + 27244*x2^2*x3^2 - 24892*x1*x3^3 + 2046*x2*x3^3 + 2516*x3^4 - 18588*x1^3*x4 + 9980*x1^2*x2*x4 - 10104*x1*x2^2*x4 + 21688*x2^3*x4 - 1315*x1^2*x3*x4 - 17824*x1*x2*x3*x4 + 14919*x2^2*x3*x4 - 568*x1*x3^2*x4 - 22509*x2*x3^2*x4 + 18494*x3^3*x4 + 25947*x1^2*x4^2 - 28652*x1*x2*x4^2 - 25547*x2^2*x4^2 + 1637*x1*x3*x4^2 - 20130*x2*x3*x4^2 + 19739*x3^2*x4^2 + 3742*x1*x4^3 + 25425*x2*x4^3 + 6342*x3*x4^3 - 3004*x4^4,
         2857*x1^4 + 8898*x1^3*x2 + 16959*x1^2*x2^2 - 28026*x1*x2^3 - 25631*x2^4 + 11030*x1^3*x3 + 29101*x1^2*x2*x3 + 30359*x1*x2^2*x3 + 27330*x2^3*x3 + 19126*x1^2*x3^2 - 26603*x1*x2*x3^2 + 2510*x2^2*x3^2 + 7575*x1*x3^3 - 25033*x2*x3^3 - 21024*x3^4 + 30501*x1^3*x4 + 23834*x1^2*x2*x4 - 1858*x1*x2^2*x4 - 10862*x2^3*x4 + 30320*x1^2*x3*x4 + 19705*x1*x2*x3*x4 + 28359*x2^2*x3*x4 + 17590*x1*x3^2*x4 + 11929*x2*x3^2*x4 + 22830*x3^3*x4 + 13501*x1^2*x4^2 - 24860*x1*x2*x4^2 + 12598*x2^2*x4^2 - 9409*x1*x3*x4^2 - 2827*x2*x3*x4^2 - 8608*x3^2*x4^2 + 30938*x1*x4^3 - 12892*x2*x4^3 + 9165*x3*x4^3 - 5202*x4^4,
         -23687*x1^4 + 32692*x1^3*x2 + 20539*x1^2*x2^2 - 27327*x1*x2^3 + 3928*x2^4 - 13018*x1^3*x3 - 13583*x1^2*x2*x3 - 30594*x1*x2^2*x3 - 12584*x2^3*x3 - 9819*x1^2*x3^2 + 14542*x1*x2*x3^2 + 30297*x2^2*x3^2 + 15188*x1*x3^3 + 28438*x2*x3^3 + 13512*x3^4 + 13327*x1^3*x4 + 14335*x1^2*x2*x4 - 15128*x1*x2^2*x4 - 21922*x2^3*x4 - 22104*x1^2*x3*x4 - 23535*x1*x2*x3*x4 - 4393*x2^2*x3*x4 - 20398*x1*x3^2*x4 + 14310*x2*x3^2*x4 - 4426*x3^3*x4 - 23087*x1^2*x4^2 - 21281*x1*x2*x4^2 + 13831*x2^2*x4^2 - 23378*x1*x3*x4^2 + 18852*x2*x3*x4^2 - 11968*x3^2*x4^2 - 31181*x1*x4^3 + 20091*x2*x4^3 - 14043*x3*x4^3 - 10677*x4^4]
    sgb = sig_groebner_basis(F)
    @test AlgebraicSolving._is_gb(sgb)
end
