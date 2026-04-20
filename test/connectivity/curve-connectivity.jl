@testset "Algorithms -> Curve Graph" begin
    R, (x,y)  = polynomial_ring(QQ, [:x,:y])

     # param2 (generic)
    f = -49303382*x^4 + 9395599*x^3*y + 67366686*x^3 - 27407214*x^2*y^2 - 71298014*x^2*y - 24150320*x^2 - 12067817*x*y^3 + 6797370*x*y^2 + 20838420*x*y - 18118613*x + 2357706*y^4 - 13736522*y^3 - 37604516*y^2 - 13868221*y + 6980802
    g = 32091920*x^4 - 97772598*x^3*y - 245256584*x^3 + 99093410*x^2*y^2 + 162335273*x^2*y + 239310556*x^2 + 28995368*x*y^3 - 59544597*x*y^2 - 20499914*x*y + 98243402*x + 22738226*y^3 + 111105506*y^2 + 86287748*y - 62439535
    P = AlgebraicSolving.RationalCurveParametrization([:z,:x,:y], Vector{Vector{ZZRingElem}}(), f, derivative(f, 2), [g])
    G = curve_graph(P)
    @test number_of_connected_components(G) == 3
end