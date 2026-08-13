@testset "Examples -> Eco" begin
    I = eco(1)
    R = first(I.gens).parent

    @test I.gens == [R(1)]

    I = eco(2)
    R, x = first(I.gens).parent, vars(first(I.gens))
    G = [x[1]*x[2] - 1,
        x[1] + 1]

    @test I.gens == G

    I = eco(5)
    R, x = first(I.gens).parent, vars(first(I.gens))
    G = [(x[1] + x[1]*x[2] + x[2]*x[3] + x[3]*x[4])*x[5] - 1,
        (x[2] + x[1]*x[3] + x[2]*x[4])*x[5] - 2,
        (x[3] + x[1]*x[4])*x[5] - 3,
        x[4]*x[5] - 4,
        x[1] + x[2] + x[3] + x[4] + 1]

    @test I.gens == G

    I = eco(4, 101)
    R, x = first(I.gens).parent, vars(first(I.gens))
    G = [(x[1] + x[1]*x[2] + x[2]*x[3])*x[4] + 100,
        (x[2] + x[1]*x[3])*x[4] + 99,
        x[3]*x[4] + 98,
        x[1] + x[2] + x[3] + 1]

    @test I.gens == G

    R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
    I = eco(R)
    G = [x*y*z + x*z - 1,
        y*z - 2,
        x + y + 1]

    @test I.gens == G

    @test_throws ArgumentError eco(0)

    R, _ = polynomial_ring(QQ, ["x"])
    @test eco(R).gens == [R(1)]
end
