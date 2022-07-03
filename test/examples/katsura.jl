@testset "Examples -> Katsura" begin
    F = katsura(3)
    R, x = first(F).parent, vars(first(F))
    G = [x[1] + 2*x[2] + 2*x[3] + 2*x[4] - 1,
         x[1]^2 + 2*x[2]^2 + 2*x[3]^2 + 2*x[4]^2 - x[1],
         2*x[1]*x[2] + 2*x[2]*x[3] + 2*x[3]*x[4] - x[2],
         x[2]^2 + 2*x[1]*x[3] + 2*x[2]*x[4] - x[3]]

    @test F == G
    
    F = katsura(4, 101)
    R, x = first(F).parent, vars(first(F))
    G = [x[1] + 2*x[2] + 2*x[3] + 2*x[4] + 2*x[5] + 100,
         x[1]^2 + 2*x[2]^2 + 2*x[3]^2 + 2*x[4]^2 + 2*x[5]^2 + 100*x[1],
         2*x[1]*x[2] + 2*x[2]*x[3] + 2*x[3]*x[4] + 2*x[4]*x[5] + 100*x[2],
         x[2]^2 + 2*x[1]*x[3] + 2*x[2]*x[4] + 2*x[3]*x[5] + 100*x[3],
         2*x[2]*x[3] + 2*x[1]*x[4] + 2*x[2]*x[5] + 100*x[4]]

    @test F == G
    
    R, (x,y,z) = PolynomialRing(QQ, ["x","y","z"])
    F = katsura(R)
    G = [x + 2*y + 2*z - 1,
         x^2 - x + 2*y^2 + 2*z^2,
         2*x*y + 2*y*z - y]

    @test F == G
end
