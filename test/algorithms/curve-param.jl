@testset "Algorithms -> Curve Parametrization" begin
    R, (x1,x2,x3) = polynomial_ring(QQ, ["x1","x2","x3"])
    I = Ideal([x1+2*x2+2*x3-1, x1^2+2*x2^2+2*x3^2-x1])

    C, (x,y) = polynomial_ring(QQ, ["x","y"])
    elim  = x^2 + 4//3*x*y - 1//3*x + y^2 - 1//3*y
    denom = 4//3*x + 2*y - 1//3

    p = [4//3*x^2 - 4//3*x*y + 2//3*x + 4//3*y - 1//3,
    -2*x^2 - 4//3*x*y + 2//3*x + 1//3*y,
    4//3*x^2 + 2*x*y - 1//3*x]

    param = curve_rational_parametrization(I)

    @test param.vars == [:x1, :x2, :x3, :_Z2, :_Z1]
    @test param.cfs_lfs == Vector{ZZRingElem}[[0, 0, 1, 0, -1], [0, 1, 0, -1, 0]]
    @test param.elim == elim
    @test param.denom == denom
    @test param.param == p
    @test I.rat_param.elim == elim
    @test I.rat_param.denom == denom
    @test I.rat_param.param == p

    R, (x1,x2,x3,x4) = polynomial_ring(QQ, ["x1","x2","x3","x4"])

    I = Ideal([x1^2-x2, x1*x3-x4, x2*x4-12, x4^3-x3^2])
    curve_rational_parametrization(I)
    @test I.rat_param.vars == Symbol[]
    @test I.rat_param.elim == -one(C)
    @test I.dim == -1

    # Automatic generic linear form
    I = Ideal([x1^2-x2, x1*x3, x2-12])

    elim    = y^2 - 12
    denom   = 2*y
    p1      = 24
    p2      = 24*y
    p3      = zero(C)
    p4      = 2*x*y

    param = curve_rational_parametrization(I)

    @test param.vars == [:x1, :x2, :x3, :x4, :_Z2, :_Z1]
    @test param.cfs_lfs == Vector{ZZRingElem}[[0, 0, 0, 1, 0, -1], [1, 0, 0, 0, -1, 0]]
    @test param.elim == elim
    @test param.denom == denom
    @test param.param[1] == p1
    @test param.param[2] == p2
    @test param.param[3] == p3
    @test param.param[4] == p4
    @test I.dim == 1

    # Prescribed generic linear form
    I = Ideal([x1^2-x2, x1*x3, x2-12])

    elim    = 2500//5041*x^2 - 100//71*x*y + 10800//5041*x + y^2 - 216//71*y - 25968//5041
    denom   = -100//71*x + 2*y - 216//71
    p1      = -1344//71
    p2      = -1200//71*x + 24*y - 2592//71
    p3      = zero(C)
    p4      = -2600//5041*x^2 + 52//71*x*y - 89616//5041*x + 1680//71*y - 20160//5041

    param = curve_rational_parametrization(I, cfs_lfs=map.(Ref(ZZ),[[-8,2,2,-1,-8,6], [8,-7,-5,8,-7,2]]))

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
	@test_throws AssertionError curve_rational_parametrization(I)

    I = Ideal([R(0)])
    @test_throws AssertionError curve_rational_parametrization(I)

    I = Ideal([R(1)])
    @test curve_rational_parametrization(I).vars == Symbol[]

    # Very non-generic variables
    A,(s,t,u) = polynomial_ring(QQ, [:s,:u,:t])
    param = curve_rational_parametrization(Ideal([u^2+t, t-s]))
    param.cfs_lfs = Vector{ZZRingElem}[[0, 1, 0, 0, -1], [1, 1, 1, -1, 0]]
    param.elim = 4*x^2 - 4*x*y + x + y^2
    param.param = QQMPolyRingElem[-4*x^2 + 2*x*y, -4*x^2 + 2*x*y, -2*x]
end
