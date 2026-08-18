@testset "Algorithms -> Semialgebraic sets" begin
    R, (x, y) = polynomial_ring(QQ, ["x", "y"])
    S = AlgebraicSolving.Intersect([
        AlgebraicSolving.BasicSemialgebraicSet(eqs=[x^2 + y^2 - 1], ineqs=[x - y], pos=[], nonneg=[]),
        AlgebraicSolving.Complement(AlgebraicSolving.BasicSemialgebraicSet(eqs=[], ineqs=[x + y - 1], pos=[], nonneg=[]))
    ])
    v = evaluate(S, [QQ(1 // 2), QQ(1 // 2)])
    @test !v
end
