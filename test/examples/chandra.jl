@testset "Examples -> Chandra" begin
    I = chandra(1)
    R, x = polynomial_ring(QQ, 1, internal_ordering=:degrevlex)
    c = QQ(51234, 100000)

    @test I.gens == [(2 - c)*x[1] - 2]

    I = chandra(3)
    R, x = polynomial_ring(QQ, 3, internal_ordering=:degrevlex)
    G = [6*x[1] - c*x[1]*(1 + QQ(1, 2)*x[1] + QQ(1, 3)*x[2]) - 6,
        6*x[2] - c*x[2]*(1 + QQ(2, 3)*x[1] + QQ(2, 4)*x[2]) - 6,
        6*x[3] - c*x[3]*(1 + QQ(3, 4)*x[1] + QQ(3, 5)*x[2]) - 6]

    @test I.gens == G

    I = chandra(3, 101)
    F = GF(101)
    R, x = polynomial_ring(F, 3, internal_ordering=:degrevlex)
    cp = F(51234) // F(100000)
    G = [6*x[1] - cp*x[1]*(1 + F(1)//F(2)*x[1] + F(1)//F(3)*x[2]) - 6,
        6*x[2] - cp*x[2]*(1 + F(2)//F(3)*x[1] + F(2)//F(4)*x[2]) - 6,
        6*x[3] - cp*x[3]*(1 + F(3)//F(4)*x[1] + F(3)//F(5)*x[2]) - 6]

    @test I.gens == G

    R, (x, y) = polynomial_ring(QQ, ["x", "y"])
    I = chandra(R, c=QQ(1, 2))
    G = [-QQ(1, 4)*x^2 + QQ(7, 2)*x - 4,
        -QQ(1, 3)*x*y + QQ(7, 2)*y - 4]

    @test I.gens == G

    @test_throws ArgumentError chandra(0)
    @test_throws ArgumentError chandra(3, 5)
    @test_throws ArgumentError chandra(3, 3)
end
