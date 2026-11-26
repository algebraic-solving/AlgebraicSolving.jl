import Nemo: fraction_field, map_coefficients, total_degree, numerator, denominator, coeff, gens, exponent_vectors, coefficients

_specialize(c, v) = evaluate(numerator(c), v) // evaluate(denominator(c), v)

@testset "Algorithms -> Chow form (univariate)" begin
    R, (y,) = polynomial_ring(QQ, [:y])
    K = fraction_field(R)
    S, (x,) = polynomial_ring(K, [:x], internal_ordering=:degrevlex)
    T, (u0, u1) = polynomial_ring(K, [:u0, :u1])
    f = (y + 1) * x^2 + y * x + y
    C = chow_form(ParametricIdeal([f]), T; reduced=false)
    @test parent(C) === T
    @test total_degree(C) == 2
    @test coeff(C, [2, 0]) == 1
    # C(1,u) = (-u)^δ f(-1/u) / lc(f)
    @test evaluate(C, [one(T), u1]) == (K(y + 1) - K(y) * u1 + K(y) * u1^2) / K(y + 1)
    @test C == chow_form(ParametricIdeal([f]), T; reduced=false, algorithm=:det)
    # the primitive Chow form lies in the same ring, with coefficients in ZZ[y]
    Cp = chow_form(ParametricIdeal([f]), T; reduced=false, primitive=true)
    @test parent(Cp) === T
    @test Cp == K(y + 1) * u0^2 - K(y) * u0 * u1 + K(y) * u1^2
    # content over ZZ[y], not only over ZZ, is removed and the sign is normalized
    @test chow_form(ParametricIdeal([-(2 * y + 2) * f]), T; reduced=false, primitive=true) == Cp
end

@testset "Algorithms -> Chow form (determinant cross-check)" begin
    R, (a, b) = polynomial_ring(QQ, [:a, :b])
    K = fraction_field(R)
    S, (x1, x2) = polynomial_ring(K, [:x1, :x2], internal_ordering=:degrevlex)
    F = [x1^2 + a * x2 - 1, x1 * x2 - b]
    I = ParametricIdeal(F)
    T, _ = polynomial_ring(K, [:v0, :v1, :v2])
    C = chow_form(I, T; reduced=false)
    @test parent(C) === T
    @test total_degree(C) == 3
    @test C == chow_form(I, T; reduced=false, algorithm=:det)
    # same result when the multiplication matrices are interpolated
    # instead of being derived from the parametric Gröbner basis
    I′ = ParametricIdeal(F; prefer_gb=false)
    @test chow_form(I′, T; reduced=false) == C
end

@testset "Algorithms -> Chow form (non-reduced)" begin
    R, (y,) = polynomial_ring(QQ, [:y])
    K = fraction_field(R)
    S, (x1, x2) = polynomial_ring(K, [:x1, :x2], internal_ordering=:degrevlex)
    T, (u0, u1, u2) = polynomial_ring(K, [:u0, :u1, :u2])
    I = ParametricIdeal([x1^2, x2 - y])
    C = chow_form(I, T; reduced=false)
    @test C == (u0 + K(y) * u2)^2
    @test C == chow_form(I, T; reduced=false, algorithm=:det)
    # a triple point and a simple point
    I = ParametricIdeal([x1^3 * (x1 - y), x2 - x1 - 1])
    C = chow_form(I, T; reduced=false)
    @test C == (u0 + u2)^3 * (u0 + K(y) * u1 + K(y + 1) * u2)
end

@testset "Algorithms -> Chow form (reduced)" begin
    R, (y,) = polynomial_ring(QQ, [:y])
    K = fraction_field(R)
    S, (x1, x2) = polynomial_ring(K, [:x1, :x2], internal_ordering=:degrevlex)
    T, (u0, u1, u2) = polynomial_ring(K, [:u0, :u1, :u2])
    # double point
    F = [x1^2, x2 - y]
    C = chow_form(ParametricIdeal(F), T)
    @test C == u0 + K(y) * u2
    @test C == chow_form(ParametricIdeal(F), T; reduced=true)
    @test C == chow_form(ParametricIdeal(F; radicalize=true), T)
    @test C == chow_form(ParametricIdeal(F), T; reduced=true, algorithm=:det)
    @test chow_form(ParametricIdeal(F), T; reduced=true, primitive=true) == C
    @test_throws ArgumentError chow_form(ParametricIdeal(F), T; algorithm=:foo)
    # the ring of the result must match the ideal
    @test_throws ArgumentError chow_form(ParametricIdeal(F), polynomial_ring(K, [:u0, :u1])[1])
    @test_throws ArgumentError chow_form(ParametricIdeal(F), polynomial_ring(fraction_field(polynomial_ring(QQ, [:z])[1]), [:u0, :u1, :u2])[1])
    # triple point and simple point
    F = [x1^3 * (x1 - y), x2 - x1 - 1]
    C = chow_form(ParametricIdeal(F), T; reduced=true)
    @test C == (u0 + u2) * (u0 + K(y) * u1 + K(y + 1) * u2)
    @test C == chow_form(ParametricIdeal(F; radicalize=true), T)
    # fat point defined by the square of the maximal ideal, of multiplicity 3
    F = [x1^2, x1 * (x2 - y), (x2 - y)^2]
    @test chow_form(ParametricIdeal(F), T; reduced=false) == (u0 + K(y) * u2)^3
    C = chow_form(ParametricIdeal(F), T; reduced=true)
    @test C == u0 + K(y) * u2
    @test C == chow_form(ParametricIdeal(F; radicalize=true), T)
    @test C == chow_form(ParametricIdeal(F), T; algorithm=:det)
    # radical ideal: same as with multiplicities
    F = [x1^2 + y * x1 * x2 - x2^2 + x1 + 2, (y + 1) * x2^2 - x1 * x2 + 3 * x1 - y]
    @test chow_form(ParametricIdeal(F), T) == chow_form(ParametricIdeal(F), T; reduced=false)
end

@testset "Algorithms -> Chow form (generic dense system)" begin
    R, (y,) = polynomial_ring(QQ, [:y])
    K = fraction_field(R)
    S, (x1, x2) = polynomial_ring(K, [:x1, :x2], internal_ordering=:degrevlex)
    F = [x1^2 + y * x1 * x2 - x2^2 + x1 + 2, (y + 1) * x2^2 - x1 * x2 + 3 * x1 - y]
    I = ParametricIdeal(F)
    T, _ = polynomial_ring(K, [:u0, :u1, :u2])
    C = chow_form(I, T; reduced=false)
    @test parent(C) === T
    δ = 4
    @test total_degree(C) == δ
    @test all(e -> sum(e) == δ, exponent_vectors(C))
    @test coeff(C, [δ, 0, 0]) == 1
    @test C == chow_form(I, T; reduced=false, algorithm=:det)
    # independent check at y = 3: C(t, c_1, c_2) is the monic eliminant of
    # t + c_1 x_1 + c_2 x_2, computed with msolve
    y₀ = [QQ(3)]
    c = [QQ(2), QQ(-3)]
    Rₜ, (z1, z2, t) = polynomial_ring(QQ, [:x1, :x2, :t])
    Fₜ = [evaluate(map_coefficients(k -> _specialize(k, y₀), f), [z1, z2]) for f in F]
    gb = AlgebraicSolving.groebner_basis(AlgebraicSolving.Ideal([Fₜ; t + c[1] * z1 + c[2] * z2]); eliminate=2)
    @test length(gb) == 1
    Cₜ = evaluate(map_coefficients(k -> _specialize(k, y₀), C), [t, Rₜ(c[1]), Rₜ(c[2])])
    @test Cₜ == gb[1] * inv(leading_coefficient(gb[1]))
    # the primitive Chow form is C multiplied by its coefficient a of u_0^δ
    Cp = chow_form(I, T; reduced=false, primitive=true)
    @test parent(Cp) === T
    a = coeff(Cp, [δ, 0, 0])
    @test Cp == a * C
    # coefficients in ZZ[y] without common factor, a with positive leading coefficient
    cs = [numerator(c) for c in coefficients(Cp)]
    @test all(c -> isone(denominator(c)), coefficients(Cp))
    @test all(p -> all(c -> isone(denominator(c)), coefficients(p)), cs)
    @test isone(reduce(gcd, cs))
    @test reduce(gcd, [numerator(c) for p in cs for c in coefficients(p)]) == 1
    @test leading_coefficient(numerator(a)) > 0
end
