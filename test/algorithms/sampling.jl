function points_per_components(ineqs::Vector{QQMPolyRingElem})::Vector{Vector{Vector{QQFieldElem}}}
    AlgebraicSolving.points_per_components(QQMPolyRingElem[], QQMPolyRingElem[], ineqs)
end

@testset "Algorithms -> Univariate sample points" begin
    R, (x,) = polynomial_ring(QQ, ["x"])
    f = [one(R)]
    @test length(points_per_components(f)) == 1
    f = [x + 1, x + 4294967295 // 4294967296, x + 4294967297 // 4294967296]
    @test length(points_per_components(f)) == 4
end

@testset "Algorithms -> Bivariate sample points" begin
    R, (x, y) = polynomial_ring(QQ, ["x", "y"])
    f = [one(R)]
    @test length(points_per_components(f)) == 1
    f = [x, x + 1]
    @test length(points_per_components(f)) == 3
    f = [y, y + 1]
    @test length(points_per_components(f)) == 3
    f = [x, x + 1, y - x, y + x - 1]
    @test length(points_per_components(f)) >= 9
end

@testset "Algorithms -> Trivariate sample points" begin
    R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
    f = [one(R)]
    @test length(points_per_components(f)) == 1
    f = [x, x + 1]
    @test length(points_per_components(f)) == 3
    f = [y, y + 1]
    @test length(points_per_components(f)) == 3
    f = [x, x + 1, y - x, y + x - 1]
    @test length(points_per_components(f)) >= 9
    f = [x, y, z]
    @test length(points_per_components(f)) == 8
end
