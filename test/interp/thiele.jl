using AlgebraicSolving.Interpolation

@testset "Interpolation -> Thiele" begin
    R, x = polynomial_ring(QQ, :x)
    for f in [zero(R), one(R), x, x * x, x // (x + 1), (x^2 + 1) // (x^4 + x^3 + 1)]
        @test thiele(R, t -> evaluate(f, t)) == f
    end
end

@testset "Interpolation -> Thiele (unlucky points)" begin
    R, x = polynomial_ring(QQ, :x)
    f = (x^2 + 1) // (x + 7)
    # The black box is degenerate on {3, 4, 6}, which Thiele's
    # interpolation hits since it evaluates at 2, 3, 4, ...
    unlucky = QQ.([3, 4, 6])
    hits = Ref(0)
    bb = function (t)
        if t in unlucky
            hits[] += 1
            throw(Interpolation.UnluckyPointError("unlucky point $t"))
        end
        evaluate(f, t)
    end
    @test thiele(R, bb) == f
    @test hits[] == length(unlucky)
    # A black box that is degenerate everywhere must still fail.
    @test_throws Interpolation.ThieleError thiele(R, t -> throw(Interpolation.UnluckyPointError("always unlucky")))
    # Genuine errors are not mistaken for unlucky points.
    @test_throws ErrorException thiele(R, t -> t == 4 ? error("genuine error") : evaluate(f, t))
end
