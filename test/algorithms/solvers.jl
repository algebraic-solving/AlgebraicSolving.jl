@testset "Algorithms -> Solvers" begin
    R, (x1,x2,x3,x4) = polynomial_ring(QQ,["x1","x2","x3","x4"], internal_ordering=:degrevlex)
    I = Ideal([x1 + 2*x2 + 2*x3 + 2*x4 - 1,
         x1^2 + 2*x2^2 + 2*x3^2 + 2*x4^2 - x1,
         2*x1*x2 + 2*x2*x3 + 2*x3*x4 - x2,
         x2^2 + 2*x1*x3 + 2*x2*x4 - x3])
    sols = Vector{QQFieldElem}[
        [4862548775//8589934592, 1281562925//8589934592, 2195068207//8589934592, -32714273694608759819673593948262790419337//174224571863520493293247799005065324265472],
        [3779635503//8589934592, 2638476131//8589934592, 908473689//8589934592, -92633843493479102248442236077889080803361//696898287454081973172991196020261297061888],
        [1, 0, 0, 0],
        [6410479475//8589934592, 2005530653//8589934592, -1585770177//8589934592, 55658687714722201275489356321691134094555601//713623846352979940529142984724747568191373312],
        [1611414365//8589934592, 673053615//8589934592, 632173751//8589934592, 708759148891639684402860468800417934359477//2787593149816327892691964784081045188247552],
        [2863311531//8589934592, 0, 0, 14518714321960041107770649917088777022123//43556142965880123323311949751266331066368]
                            ]
    rat_sols = Vector{QQFieldElem}[[49, 0, 0, 0], [49//3, 0, 0, 1//3]]

    @test sols == real_solutions(I)
    @test rat_sols == rational_solutions(I)
    @test I.real_sols == real_solutions(I)
    
    C, x = polynomial_ring(QQ, "x")
    elim  = 128304*x^8 - 93312*x^7 + 15552*x^6 + 3144*x^5 - 1120*x^4 + 36*x^3 + 15*x^2 - x
    denom = 1026432*x^7 - 653184*x^6 + 93312*x^5 + 15720*x^4 - 4480*x^3 + 108*x^2 + 30*x - 1
    p1    = -3872448*x^7 + 2607552*x^6 - 408528*x^5 - 63088*x^4 + 20224*x^3 - 540*x^2 - 172*x + 7
    p2    = -303264*x^7 + 314928*x^6 - 113544*x^5 + 9840*x^4 + 3000*x^3 - 564*x^2 + 12*x
    p3    = -699840*x^7 + 449712*x^6 - 74808*x^5 - 1956*x^4 + 1308*x^3 - 174*x^2 + 18*x
    p1   *= -7
    p2   *= -7
    p3   *= -7

    param = rational_parametrization(I)

    @test param.vars == [:x1, :x2, :x3, :x4]
    @test param.lf_cfs == ZZRingElem[]
    @test param.elim == elim
    @test param.denom == denom
    @test param.param[1] == p1
    @test param.param[2] == p2
    @test param.param[3] == p3
    @test I.rat_param.elim == elim
    @test I.rat_param.denom == denom
    @test I.rat_param.param[1] == p1
    @test I.rat_param.param[2] == p2
    @test I.rat_param.param[3] == p3

    I = Ideal([x1^2-x2, x1*x3-x4, x2*x4-12, x4^3-x3^2])
    real_solutions(I)
    @test I.rat_param.vars == Symbol[]

    I = Ideal([x1^2-x2, x1*x3, x2-12])
	@test_throws ErrorException real_solutions(I)
	@test_throws ErrorException rational_solutions(I)
end
