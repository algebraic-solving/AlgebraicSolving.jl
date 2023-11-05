@testset "Algorithms -> Normal forms" begin
    R, (x,y) = polynomial_ring(GF(101),["x","y"])
    I = Ideal([y*x+17-y, x+13*y])
    f = x^4+y
    @test normal_form(f, I) == 100*y+77
    F = [x+13*y, x*y-4]
    @test normal_form(F, I) == [0, y]
end
