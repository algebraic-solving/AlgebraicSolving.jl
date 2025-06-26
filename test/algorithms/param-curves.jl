@testset "Algorithms -> Curve Parametrization" begin
    R, (x1,x2,x3) = polynomial_ring(QQ, ["x1","x2","x3"])
    I = Ideal([x1+2*x2+2*x3-1, x1^2+2*x2^2+2*x3^2-x1])

    C, (x,y) = polynomial_ring(QQ, ["x","y"])
    elim  = x^2 + 4//3*x*y - 1//3*x + y^2 - 1//3*y
    denom = 4//3*x + 2*y - 1//3
    p     = 4//3*x^2 - 4//3*x*y + 2//3*x + 4//3*y - 1//3

    param = rational_curve_parametrization(I)

    @test param.vars == [:x1, :x2, :x3]
    @test param.cfs_lfs == Vector{ZZRingElem}[]
    @test param.elim == elim
    @test param.denom == denom
    @test param.param[1] == p1
    @test I.rat_param.elim == elim
    @test I.rat_param.denom == denom
    @test I.rat_param.param[1] == p1

    R, (x1,x2,x3,x4) = polynomial_ring(QQ, ["x1","x2","x3","x4"])

    I = Ideal([x1^2-x2, x1*x3-x4, x2*x4-12, x4^3-x3^2])
    rational_curve_parametrization(I)
    @test I.rat_param.vars == Symbol[]
    @test I.rat_param.elim == -one(C)
    @test I.dim == -1

    I = Ideal([x1^2-x2, x1*x3, x2-12])
	@test_throws AssertionException rational_curve_parametrization(I)
    @test_throws ErrorException rational_curve_parametrization(Ideal([x1^2-x2, x1*x3, x2-12]), check_gen=false)

    C, (x,y,_z1,_z2) = polynomial_ring(QQ, ["x","y","_Z2","_Z1"])
    elim    = 5041//15129*x^2 + 142//123*x*y - 1704//1681*x + y^2 - 72//41*y - 8656//5043
    denom   = 142//123*x + 2*y - 72//41
    p1      = -448//41
    p2      = 568//41*x + 24*y - 864//41
    p3      = zero(C)
    p4      = 13916//5043*x^2 + 196//41*x*y + 8848//1681*x + 672//41*y - 27776//1681

    param = rational_curve_parametrization(I, cfs_lfs=map.(Ref(ZZ),[[-8,2,2,-1,-8], [8,-7,-5,8,-7]]))

    @test param.vars == [:x1, :x2, :x3, :x4, :_Z2, :_Z1]
    @test param.cfs_lfs == Vector{ZZRingElem}[[-8, 2, 2, -1, -8], [8, -7, -5, 8, -7]]
    @test param.elim == elim
    @test param.denom == denom
    @test param.param[1] == p1
    @test param.param[2] == p2
    @test param.param[3] == p3
    @test param.param[4] == p4
    @test I.dim = 1

    I = Ideal([x1^2-x2, x1*x3])
	@test_throws AssertionException rational_curve_parametrization(I)
    @test_throws AssertionException rational_curve_parametrization(I, use_lfs=true)
    @test_throws AssertionException rational_curve_parametrization(I, check_gen=false)

    I = Ideal([R(0)])
    @test_throws AssertionException rational_curve_parametrization(I)

    I = Ideal([R(1)])
    @test rational_curve_parametrization(I).vars == Symbol[]
end