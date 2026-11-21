using AlgebraicSolving.Interpolation

@testset "Interpolation -> Newton" begin
    R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
    n = 5
    xs = [ZZ(i) for i in 1:n]
    for f in [zero(R), one(R), x, y, z, x + z, y * z, x * y * (x + y + z)]
        @test newton(R, t -> evaluate(f, t), [n - 1, n - 1, n - 1]) == f
    end
end

@testset "Interpolation -> Newton (total degree)" begin
    I = AlgebraicSolving.Interpolation
    # the simplex lattice has binomial(n + d, n) points, graded by total degree
    for n in 0:4, d in 0:5
        idx = I._simplex_indices(n, d)
        @test length(idx) == binomial(n + d, n)
        @test allunique(idx)
        @test issorted(sum.(idx))
        @test all(a -> length(a) == n && all(>=(0), a) && sum(a) <= d, idx)
    end

    R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
    for f in [zero(R), one(R), x, y, z, x + z, y * z, x * y * (x + y + z), x^4 - 3 * y^2 * z^2 // 7 + z - 1]
        @test newton(R, t -> evaluate(f, t), 4) == f
        @test newton(R, t -> evaluate(f, t), 6) == f
        @test newton(R, t -> evaluate(f, t), 4) == newton(R, t -> evaluate(f, t), [4, 4, 4])
    end
    @test newton(R, t -> evaluate(x * y * z, t), 3; nr_thrds=4) == x * y * z

    # no variables / one variable
    R0, _ = polynomial_ring(QQ, Symbol[])
    @test newton(R0, t -> QQ(5), 3) == R0(5)
    R1, (u,) = polynomial_ring(QQ, [:u])
    @test newton(R1, t -> evaluate(u^5 - u // 3, t), 5) == u^5 - u // 3

    # the black box is never evaluated at zeros of non_vanishing_poly
    f = x^3 * y - 7 * z^2 + x * y * z + 1
    for p in [x - y, (x - 1) * (y - 1), x - 2 * y + z, (x - y) * (y - z) * (x - z), x * y * z - 6]
        bb = t -> (@assert !iszero(evaluate(p, t)); evaluate(f, t))
        @test newton(R, bb, 4; non_vanishing_poly=p) == f
    end
end

@testset "Interpolation -> Newton (node choice)" begin
    I = AlgebraicSolving.Interpolation
    R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
    d = 5
    idx = I._simplex_indices(3, d)
    # default nodes when nothing vanishes on the default lattice
    @test I._choose_nodes(one(R), d) == [ZZ.(1:d+1) for _ in 1:3]
    @test I._choose_nodes(x + y + z + 100, d) == [ZZ.(1:d+1) for _ in 1:3]
    # p never vanishes on the lattice, nodes are distinct positive integers on each axis
    for p in [x - y, (x - 1) * (y - 1), x - 2 * y + z, (x - y) * (y - z) * (x - z), x * y * z - 6,
              prod(x - c for c in 1:6), (x - y)^2 * (y - 3), (x + y - 5) * (x * y - 6) * (z - x),
              (x - 1) * (y - 1) * (z - 1) * (x + y + z - 9)]
        nodes = I._choose_nodes(p, d)
        @test all(nj -> length(nj) == d + 1 && allunique(nj) && all(>(0), nj), nodes)
        @test all(a -> !iszero(evaluate(p, QQ.(I._simplex_point(nodes, a)))), idx)
    end
    # the x-nodes must skip the roots 1..6 of the first factor
    @test I._choose_nodes(prod(x - c for c in 1:6), d)[1] == ZZ.(7:12)
end
