@testset "Algorithms -> Normal forms" begin
    R, (x,y) = polynomial_ring(GF(101),["x","y"])
    I = Ideal([y*x+17-y, x+13*y])
    f = x^4+y
    @test normal_form(f, I) == 100*y+77
    F = [x+13*y, x*y-4]
    @test normal_form(F, I) == [0, y + 80]
    R, (x,y) = polynomial_ring(GF(536870923),["x","y"])
    I = Ideal([y*x+17-y, x+13*y])
    F = [x+13*y, x*y-4]
    @test normal_form(F, I) == [0, y + 536870902]
    R, (x,y) = polynomial_ring(GF(65521),["x","y"])
    I = Ideal([y*x+17-y, x+13*y])
    F = [x+13*y, x*y+16]
    @test normal_form(F, I) == [0, y + 65520]
    I = Ideal([R(0)])
    @test normal_form(F, I) == F
end
