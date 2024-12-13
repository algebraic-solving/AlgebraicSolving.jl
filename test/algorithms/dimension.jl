@testset "Algorithms -> Dimension" begin
    R, (x,y) = polynomial_ring(QQ,["x","y"])
    I =  Ideal([x^2,x*y])
    @test isone(dimension(I))

    R, (x,y,z) = polynomial_ring(GF(101),["x","y","z"])
    I = Ideal([x+2*y+2*z-1, x^2+2*y^2+2*z^2-x, 2*x*y+2*y*z-y])
    @test iszero(dimension(I))

    I = Ideal([R(1)])
    @test dimension(I) == -1
end
