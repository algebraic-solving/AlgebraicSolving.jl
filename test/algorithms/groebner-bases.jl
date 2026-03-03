@testset "Algorithms -> Gröbner bases" begin
    R, (x,y,z) = polynomial_ring(GF(101),["x","y","z"], internal_ordering=:degrevlex)
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

    I = Ideal([x+2*y+2*z-1, x^2+2*y^2+2*z^2-x, 2*x*y+2*y*z-y])
    G = groebner_basis(I, eliminate=2, intersect=false)
    H = MPolyRingElem[
        z^4 + 38*z^3 + 95*z^2 + 95*z
        30*z^3 + 32*z^2 + y + 87*z
        41*z^3 + 37*z^2 + x + 30*z + 100
        ]
    @test G == H

    @test_throws ErrorException eliminate(I,0)
    L = eliminate(I,2)
    @test L == H
    @test I.gb[2] == H

    I = Ideal([R(0)])
    G = groebner_basis(I)
    @test G == [R(0)]

    R, (x1, x2) = polynomial_ring(QQ, ["x1", "x2"])
    I = Ideal([3*x1^2 + ZZRingElem(2)^100, 2*x1*x2 + 5*x1 + ZZRingElem(2)^100])
    G = groebner_basis(I)
    H = MPolyRingElem[
        3*x1 - 2*x2 - 5
        4*x2^2 + 20*x2 + 3802951800684688204490109616153
        ]
    @test G == H
    J = Ideal([3*x1^2 + ZZRingElem(2)^100, 2*x1*x2 + 5*x1 + ZZRingElem(2)^100])
    G = groebner_basis(J, normalize=true)
    H = MPolyRingElem[
        x1 - 2//3*x2 - 5//3
        x2^2 + 5*x2 + 3802951800684688204490109616153//4
        ]
    @test G == H
    R, (x,y,z) = polynomial_ring(QQ,["x","y","z"], internal_ordering=:degrevlex)
    I = Ideal([x+2*y+2*z-1, x^2+2*y^2+2*z^2-x, 2*x*y+2*y*z-y])
    G = eliminate(I, 2)
    H = MPolyRingElem[
        84*z^4 - 40*z^3 + z^2 + z
    ]
    @test G == H
    R, (a,b,c,d,x,y,z,w) = polynomial_ring(QQ, ["a", "b", "c", "d", "x", "y", "z", "w"])
    I = Ideal([x - a*c, y - a*c*d, z - a*c^2 - b, w - a*c^2*d - b*d]);
    G = eliminate(I, 4)
    H = MPolyRingElem[
        -x*w + y*z
    ]
    @test G == H

		# issue 113
		R, (u1,u2,u3,u4,u5,u6,u7,u8,u9,l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,l15,l16,l17,l18,x1,x2) = polynomial_ring(QQ, [:u1, :u2, :u3, :u4, :u5, :u6, :u7, :u8, :u9, :l1, :l2, :l3, :l4, :l5, :l6, :l7, :l8, :l9, :l10, :l11, :l12, :l13, :l14, :l15, :l16, :l17, :l18, :x1, :x2])
		I = Ideal([
           -95 * u1 * x1 + 53 * u4 * x1 + 63 * u1 * x2 + 89 * u4 * x2 + u1,
           -95 * u2 * x1 + 53 * u5 * x1 + 63 * u2 * x2 + 89 * u5 * x2 + u2,
           -95 * u3 * x1 + 53 * u6 * x1 + 63 * u3 * x2 + 89 * u6 * x2 + u3,
           53 * u1 * x1 + 95 * u4 * x1 + 89 * u1 * x2 - 63 * u4 * x2 + u4,
           53 * u2 * x1 + 95 * u5 * x1 + 89 * u2 * x2 - 63 * u5 * x2 + u5,
           53 * u3 * x1 + 95 * u6 * x1 + 89 * u3 * x2 - 63 * u6 * x2 + u6,
           -95 * u7 * x1 + 63 * u7 * x2,
           -95 * u8 * x1 + 63 * u8 * x2,
           -95 * u9 * x1 + 63 * u9 * x2,
           -23 * u1 + 9 * u4 + 79 * u7 + 33,
           -23 * u2 + 9 * u5 + 79 * u8 - 22,
           -23 * u3 + 9 * u6 + 79 * u9 - 19,
           21 * u1 - 80 * u4 - 76 * u7 + 57,
           21 * u2 - 80 * u5 - 76 * u8 + 97,
           21 * u3 - 80 * u6 - 76 * u9 + 78,
           46 * u1 - 50 * u4 + 28 * u7 - 7,
           46 * u2 - 50 * u5 + 28 * u8 - 32,
           46 * u3 - 50 * u6 + 28 * u9 + 29,
           -95 * u1 * l1 + 53 * u4 * l1 - 95 * u2 * l2 + 53 * u5 * l2 - 95 * u3 * l3 + 53 * u6 * l3 + 53 * u1 * l4 + 95 * u4 * l4 + 53 * u2 * l5 + 95 * u5 * l5 + 53 * u3 * l6 + 95 * u6 * l6 - 95 * u7 * l7 - 95 * u8 * l8 - 95 * u9 * l9 - 1,
           63 * u1 * l1 + 89 * u4 * l1 + 63 * u2 * l2 + 89 * u5 * l2 + 63 * u3 * l3 + 89 * u6 * l3 + 89 * u1 * l4 - 63 * u4 * l4 + 89 * u2 * l5 - 63 * u5 * l5 + 89 * u3 * l6 - 63 * u6 * l6 + 63 * u7 * l7 + 63 * u8 * l8 + 63 * u9 * l9,
           -95 * l1 * x1 + 53 * l4 * x1 + 63 * l1 * x2 + 89 * l4 * x2 + l1 - 23 * l10 + 21 * l13 + 46 * l16,
           -95 * l2 * x1 + 53 * l5 * x1 + 63 * l2 * x2 + 89 * l5 * x2 + l2 - 23 * l11 + 21 * l14 + 46 * l17,
           -95 * l3 * x1 + 53 * l6 * x1 + 63 * l3 * x2 + 89 * l6 * x2 + l3 - 23 * l12 + 21 * l15 + 46 * l18,
           53 * l1 * x1 + 95 * l4 * x1 + 89 * l1 * x2 - 63 * l4 * x2 + l4 + 9 * l10 - 80 * l13 - 50 * l16,
           53 * l2 * x1 + 95 * l5 * x1 + 89 * l2 * x2 - 63 * l5 * x2 + l5 + 9 * l11 - 80 * l14 - 50 * l17,
           53 * l3 * x1 + 95 * l6 * x1 + 89 * l3 * x2 - 63 * l6 * x2 + l6 + 9 * l12 - 80 * l15 - 50 * l18,
           -95 * l7 * x1 + 63 * l7 * x2 + 79 * l10 - 76 * l13 + 28 * l16,
           -95 * l8 * x1 + 63 * l8 * x2 + 79 * l11 - 76 * l14 + 28 * l17,
           -95 * l9 * x1 + 63 * l9 * x2 + 79 * l12 - 76 * l15 + 28 * l18
       ])
		G = groebner_basis(I, eliminate=27)
    H = MPolyRingElem[
        R(1)
    ]
    @test G == H
	
end

@testset "Algorithms -> Sig Gröbner bases" begin
    R, (x,y,z) = polynomial_ring(QQ,["x","y","z"], internal_ordering=:degrevlex)
    F = [x^2+1-3, x*y-z, x*z^2-3*y^2]
    #= not a finite field =#
    @test_throws ErrorException sig_groebner_basis(F)
    R, (x,y,z) = polynomial_ring(GF(17),["x","y","z"], internal_ordering=:degrevlex)
    F = [x^2+1-3, x*y-z, x*z^2-3*y^2]
    #= not homogeneous =#
    @test_throws ErrorException sig_groebner_basis(F, mod_ord = :DPOT)

    #= GB test 1 =#
    sgb = sig_groebner_basis(F)
    @test AlgebraicSolving._is_gb(sgb)
    Fhom = homogenize(F)
    sgb = sig_groebner_basis(Fhom, mod_ord = :DPOT)
    @test AlgebraicSolving._is_gb(sgb)

    #= GB test 2 =#
    R, (x,y,z,w) = polynomial_ring(GF(65521),["x","y","z","w"], internal_ordering=:degrevlex)
    F = cyclic(R).gens
    sgb = sig_groebner_basis(F)
    @test AlgebraicSolving._is_gb(sgb)
    Fhom = homogenize(F)
    sgb = sig_groebner_basis(Fhom, mod_ord = :DPOT)
    @test AlgebraicSolving._is_gb(sgb)

    #= GB test 3 =#
    F = katsura(R).gens
    sgb = sig_groebner_basis(F)
    @test AlgebraicSolving._is_gb(sgb)
    Fhom = homogenize(F)
    sgb = sig_groebner_basis(Fhom, mod_ord = :DPOT)
    @test AlgebraicSolving._is_gb(sgb)

    #= GB test 4 (pivot setting bug) =#
    R, (x1, x2, x3, x4) = polynomial_ring(GF(65521), ["x1", "x2", "x3", "x4"], internal_ordering=:degrevlex)
    F = [11523*x1^4 + 30378*x1^3*x2 + 30154*x1^2*x2^2 + 10157*x1*x2^3 - 28136*x2^4 - 4771*x1^3*x3 - 21056*x1^2*x2*x3 + 15696*x1*x2^2*x3 - 16144*x2^3*x3 - 1553*x1^2*x3^2 - 30379*x1*x2*x3^2 - 12735*x2^2*x3^2 + 18058*x1*x3^3 + 24670*x2*x3^3 - 16379*x3^4 + 24196*x1^3*x4 - 19411*x1^2*x2*x4 + 17610*x1*x2^2*x4 - 5715*x2^3*x4 - 21186*x1^2*x3*x4 - 22865*x1*x2*x3*x4 - 1939*x2^2*x3*x4 - 5685*x1*x3^2*x4 + 8508*x2*x3^2*x4 + 21819*x3^3*x4 - 24868*x1^2*x4^2 - 18233*x1*x2*x4^2 - 14116*x2^2*x4^2 + 28291*x1*x3*x4^2 - 9068*x2*x3*x4^2 - 15138*x3^2*x4^2 + 8921*x1*x4^3 - 18808*x2*x4^3 - 3005*x3*x4^3 + 7368*x4^4,
         31703*x1^4 + 23616*x1^3*x2 + 20696*x1^2*x2^2 - 7125*x1*x2^3 + 15334*x2^4 + 26619*x1^3*x3 + 2173*x1^2*x2*x3 - 31312*x1*x2^2*x3 - 31386*x2^3*x3 - 25244*x1^2*x3^2 - 28729*x1*x2*x3^2 + 27244*x2^2*x3^2 - 24892*x1*x3^3 + 2046*x2*x3^3 + 2516*x3^4 - 18588*x1^3*x4 + 9980*x1^2*x2*x4 - 10104*x1*x2^2*x4 + 21688*x2^3*x4 - 1315*x1^2*x3*x4 - 17824*x1*x2*x3*x4 + 14919*x2^2*x3*x4 - 568*x1*x3^2*x4 - 22509*x2*x3^2*x4 + 18494*x3^3*x4 + 25947*x1^2*x4^2 - 28652*x1*x2*x4^2 - 25547*x2^2*x4^2 + 1637*x1*x3*x4^2 - 20130*x2*x3*x4^2 + 19739*x3^2*x4^2 + 3742*x1*x4^3 + 25425*x2*x4^3 + 6342*x3*x4^3 - 3004*x4^4,
         2857*x1^4 + 8898*x1^3*x2 + 16959*x1^2*x2^2 - 28026*x1*x2^3 - 25631*x2^4 + 11030*x1^3*x3 + 29101*x1^2*x2*x3 + 30359*x1*x2^2*x3 + 27330*x2^3*x3 + 19126*x1^2*x3^2 - 26603*x1*x2*x3^2 + 2510*x2^2*x3^2 + 7575*x1*x3^3 - 25033*x2*x3^3 - 21024*x3^4 + 30501*x1^3*x4 + 23834*x1^2*x2*x4 - 1858*x1*x2^2*x4 - 10862*x2^3*x4 + 30320*x1^2*x3*x4 + 19705*x1*x2*x3*x4 + 28359*x2^2*x3*x4 + 17590*x1*x3^2*x4 + 11929*x2*x3^2*x4 + 22830*x3^3*x4 + 13501*x1^2*x4^2 - 24860*x1*x2*x4^2 + 12598*x2^2*x4^2 - 9409*x1*x3*x4^2 - 2827*x2*x3*x4^2 - 8608*x3^2*x4^2 + 30938*x1*x4^3 - 12892*x2*x4^3 + 9165*x3*x4^3 - 5202*x4^4,
         -23687*x1^4 + 32692*x1^3*x2 + 20539*x1^2*x2^2 - 27327*x1*x2^3 + 3928*x2^4 - 13018*x1^3*x3 - 13583*x1^2*x2*x3 - 30594*x1*x2^2*x3 - 12584*x2^3*x3 - 9819*x1^2*x3^2 + 14542*x1*x2*x3^2 + 30297*x2^2*x3^2 + 15188*x1*x3^3 + 28438*x2*x3^3 + 13512*x3^4 + 13327*x1^3*x4 + 14335*x1^2*x2*x4 - 15128*x1*x2^2*x4 - 21922*x2^3*x4 - 22104*x1^2*x3*x4 - 23535*x1*x2*x3*x4 - 4393*x2^2*x3*x4 - 20398*x1*x3^2*x4 + 14310*x2*x3^2*x4 - 4426*x3^3*x4 - 23087*x1^2*x4^2 - 21281*x1*x2*x4^2 + 13831*x2^2*x4^2 - 23378*x1*x3*x4^2 + 18852*x2*x3*x4^2 - 11968*x3^2*x4^2 - 31181*x1*x4^3 + 20091*x2*x4^3 - 14043*x3*x4^3 - 10677*x4^4]
    sgb = sig_groebner_basis(F)
    @test AlgebraicSolving._is_gb(sgb)
end
