@testset "Algorithms -> Curve Parametrization" begin
    # ----------------------------------------------------
    # 1. Standard intersecting quadrics
    # ----------------------------------------------------
    R, (x1, x2, x3) = polynomial_ring(QQ, ["x1", "x2", "x3"])
    I = AlgebraicSolving.Ideal([x1 + 2*x2 + 2*x3 - 1, x1^2 + 2*x2^2 + 2*x3^2 - x1])

    C, (x, y) = polynomial_ring(QQ, ["x", "y"])
    elim  = x^2 + 4//3*x*y - 1//3*x + y^2 - 1//3*y
    denom = 4//3*x + 2*y - 1//3
    p = [
        4//3*x^2 - 4//3*x*y + 2//3*x + 4//3*y - 1//3,
        -2*x^2 - 4//3*x*y + 2//3*x + 1//3*y,
        4//3*x^2 + 2*x*y - 1//3*x
    ]

    param = curve_rational_parametrization(I)

    @test param.vars == [:x1, :x2, :x3, :_Z2, :_Z1]
    @test param.cfs_lfs == Vector{ZZRingElem}[[0, 0, 1, 0, -1], [0, 1, 0, -1, 0]]
    @test param.elim == elim
    @test param.denom == denom
    @test param.param == p
    @test I.rat_param.elim == elim

    # ----------------------------------------------------
    # 2. Empty set (Dimension -1)
    # ----------------------------------------------------
    R, (x1, x2, x3, x4) = polynomial_ring(QQ, ["x1", "x2", "x3", "x4"])
    I_empty = AlgebraicSolving.Ideal([x1^2 - x2, x1*x3 - x4, x2*x4 - 12, x4^3 - x3^2])
    curve_rational_parametrization(I_empty)
    @test I_empty.rat_param.vars == Symbol[]
    @test I_empty.rat_param.elim == -one(C)
    @test I_empty.dim == -1

    # ----------------------------------------------------
    # 3. Automatic generic linear forms (Dimension 1)
    # ----------------------------------------------------
    I = AlgebraicSolving.Ideal([x1^2 - x2, x1*x3, x2 - 12])
    elim    = y^2 - 12
    denom   = 2*y
    p1, p2, p3, p4 = C(24), 24*y, zero(C), 2*x*y

    param = curve_rational_parametrization(I)

    @test param.vars == [:x1, :x2, :x3, :x4, :_Z2, :_Z1]
    @test param.cfs_lfs == Vector{ZZRingElem}[[0, 0, 0, 1, 0, -1], [1, 0, 0, 0, -1, 0]]
    @test param.elim == elim
    @test param.denom == denom
    @test param.param == [p1, p2, p3, p4]
    @test I.dim == 1

    # ----------------------------------------------------
    # 4. Prescribed valid generic linear forms
    # ----------------------------------------------------
    elim    = 2500//5041*x^2 - 100//71*x*y + 10800//5041*x + y^2 - 216//71*y - 25968//5041
    denom   = -100//71*x + 2*y - 216//71
    p1      = C(-1344//71)
    p2      = -1200//71*x + 24*y - 2592//71
    p3      = zero(C)
    p4      = -2600//5041*x^2 + 52//71*x*y - 89616//5041*x + 1680//71*y - 20160//5041

    param = curve_rational_parametrization(I, cfs_lfs=[[-8,2,2,-1,-8,6], [8,-7,-5,8,-7,2]])

    @test param.elim == elim
    @test param.denom == denom
    @test param.param == [p1, p2, p3, p4]

    # ----------------------------------------------------
    # 5. Guard Assertions & Custom Error Messages
    # ----------------------------------------------------
    # Improper sizing guards
    @test_throws AssertionError curve_rational_parametrization(I, cfs_lfs=[[1,2,3,4], [1,2,3,4,5]])
    @test_throws AssertionError curve_rational_parametrization(I, cfs_lfs=[[1,2,3,4]])

    # Dimension checking guards
    I_high_dim = AlgebraicSolving.Ideal([x1^2 - x2, x1*x3])
    @test_throws AssertionError curve_rational_parametrization(I_high_dim)
    @test_throws AssertionError curve_rational_parametrization(AlgebraicSolving.Ideal([R(0)]))

    # Non-generic linear forms
    @test_throws ["failed", "number 2"] curve_rational_parametrization(I, cfs_lfs=[[1,2,3,4], [1,2,3,4]]) # identical (non-generic)
    @test_throws ["failed", "number 1"] curve_rational_parametrization(I, cfs_lfs=[[0,0,1,0], [0,0,0,1]]) # first is not generic
    @test_throws ["failed", "number 2"] curve_rational_parametrization(I, cfs_lfs=[[0,0,0,1], [0,0,1,0]]) # second is not generic

    # Minimal edge cases
    I_one = AlgebraicSolving.Ideal([R(1)])
    @test curve_rational_parametrization(I_one).vars == Symbol[]

    # ----------------------------------------------------
    # 6. Very non-generic variables (Fixed assignment bug)
    # ----------------------------------------------------
    I_nongen = AlgebraicSolving.Ideal([x1^2 - x2, x1*x3, x4])
    param = curve_rational_parametrization(I_nongen)

    @test param.cfs_lfs == Vector{ZZRingElem}[[1, 1, 1, 1, 0, -1], [2, 1, 1, 1, -1, 0]]
    @test param.elim == -x^3 + 3*x^2*y + 2*x^2 - 3*x*y^2 - 3*x*y + y^3 + y^2
    @test param.param == QQMPolyRingElem[
        -3*x^2 + 4*x*y - y^2,
        2*x^3 - 4*x^2*y + 2*x^2 + 2*x*y^2 - 3*x*y + y^2,
        x^3 - 2*x^2*y - 2*x^2 + x*y^2 + x*y,
        C(0)
    ]

    # Non radical curve
    @test_throws ["bad specializations", "radicality"] curve_rational_parametrization(AlgebraicSolving.Ideal([x1^2, x2, x3]))
end
