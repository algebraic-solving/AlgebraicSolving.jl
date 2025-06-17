@testset "Algorithms -> Hilbert" begin
    R, (x,y) = polynomial_ring(QQ,["x","y"])
    I =  Ideal([x^2,x*y])
    A, t = polynomial_ring(ZZ, :t)
    B, s = polynomial_ring(QQ, :s)
    HS = (t^2 - t - 1)//(t - 1)
    HP = (s + 1, 2)

    @test HS == hilbert_series(I)
    @test HP == hilbert_polynomial(I)
    @test isone(hilbert_dimension(I))
    @test isone(hilbert_degree(I))

    R, (x,y,z) = polynomial_ring(GF(101),["x","y","z"])
    I = Ideal([x+2*y+2*z-1, x^2+2*y^2+2*z^2-x, 2*x*y+2*y*z-y])
    HS = t^2 + 2*t + 1
    HP = (s^2 + 2*s + 1, 0)

    @test HS == hilbert_series(I)
    @test HP == hilbert_polynomial(I)
    @test iszero(hilbert_dimension(I))
    @test 4 == hilbert_degree(I)

    I = Ideal([R(0)])
    HS = 1//(1-t)^(nvars(R))
    HP = (1//6*s^3 + s^2 + 11//6*s + 1, 0)

    @test HS == hilbert_series(I)
    @test HP == hilbert_polynomial(I)
    @test nvars(R) == hilbert_dimension(I)
    @test isone(hilbert_degree(I))

    I = Ideal([R(1)])
    @test iszero(hilbert_series(I))
    @test all(iszero, hilbert_polynomial(I))
    @test hilbert_dimension(I) == -1
    @test iszero(hilbert_degree(I))
end
