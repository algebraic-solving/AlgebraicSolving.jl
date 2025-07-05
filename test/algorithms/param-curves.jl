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
    @test param.param[1] == p
    @test I.rat_param.elim == elim
    @test I.rat_param.denom == denom
    @test I.rat_param.param[1] == p

    R, (x1,x2,x3,x4) = polynomial_ring(QQ, ["x1","x2","x3","x4"])

    I = Ideal([x1^2-x2, x1*x3-x4, x2*x4-12, x4^3-x3^2])
    rational_curve_parametrization(I)
    @test I.rat_param.vars == Symbol[]
    @test I.rat_param.elim == -one(C)
    @test I.dim == -1

    I = Ideal([x1^2-x2, x1*x3, x2-12])
	@test_throws AssertionError rational_curve_parametrization(I)
    @test_throws ErrorException rational_curve_parametrization(Ideal([x1^2-x2, x1*x3, x2-12]), check_gen=false)

    elim    = 5041//2500*x^2 - 71//25*x*y - 3834//625*x + y^2 + 108//25*y - 6492//625
    denom   = -71//25*x + 2*y + 108//25
    p1      = 672//25
    p2      = -852//25*x + 24*y + 1296//25
    p3      = zero(C)
    p4      = -923//625*x^2 + 26//25*x*y - 18192//625*x + 552//25*y + 8304//625

    param = rational_curve_parametrization(I, cfs_lfs=map.(Ref(ZZ),[[-8,2,2,-1,-8,6], [8,-7,-5,8,-7,2]]))

    @test param.vars == [:x1, :x2, :x3, :x4, :_Z2, :_Z1]
    @test param.cfs_lfs == Vector{ZZRingElem}[[-8, 2, 2, -1, -8, 6], [8, -7, -5, 8, -7, 2]]
    @test param.elim == elim
    @test param.denom == denom
    @test param.param[1] == p1
    @test param.param[2] == p2
    @test param.param[3] == p3
    @test param.param[4] == p4
    @test I.dim == 1

    I = Ideal([x1^2-x2, x1*x3])
	@test_throws AssertionError rational_curve_parametrization(I)
    @test_throws AssertionError rational_curve_parametrization(I, use_lfs=true)
    @test_throws AssertionError rational_curve_parametrization(I, check_gen=false)

    I = Ideal([R(0)])
    @test_throws AssertionError rational_curve_parametrization(I)

    I = Ideal([R(1)])
    @test rational_curve_parametrization(I).vars == Symbol[]
end