Var = AlgebraicSolving.sign_variations
Sign = AlgebraicSolving.signature
Mat = AlgebraicSolving.matrix_of_signs
Ada = AlgebraicSolving.adapted_family
∧ = AlgebraicSolving.extend_signs
× = AlgebraicSolving.extend_indices
⊗′ = AlgebraicSolving.modified_tensor_product

@testset "Algorithms -> Sign variations" begin
    @test Var([]) == 0
    @test Var([1, -1, 2, 0, 0, 3, 4, -5, -2, 0, 3]) == 4
end

@testset "Algorithms -> Signature of quadratic forms" begin
    @test Sign(matrix(QQ, [1;;])) == 1
    @test Sign(matrix(QQ, [1 2; 2 3])) == 0
    @test Sign(matrix(QQ, [2 0; 0 0])) == 1
end

@testset "Algorithms -> Matrix of signs" begin
    Σ = [Vector{Int}()] ∧ [0, 1, -1]
    A = [Vector{Int}()] × [0, 1, 2]
    @test Σ == [[0], [1], [-1]]
    @test A == [[0], [1], [2]]
    M = M′ = Mat(A, Σ)
    @test M == [
        1 1 1
        0 1 -1
        0 1 1
    ]
    @test M⊗′M′ == [
        1 1 1 1 1 1 1 1 1;
        0 0 0 1 1 1 -1 -1 -1;
        0 0 0 1 1 1 1 1 1;
        0 1 -1 0 1 -1 0 1 -1;
        0 0 0 0 1 -1 0 -1 1;
        0 0 0 0 1 -1 0 1 -1;
        0 1 1 0 1 1 0 1 1;
        0 0 0 0 1 1 0 -1 -1;
        0 0 0 0 1 1 0 1 1
    ]
    Σ = Σ ∧ [0, 1, -1]
    A = A × [0, 1, 2]
    @test Mat(A, Σ) == M⊗′M′
end

# Example 10.74
@testset "Algorithms -> Adapted family 1" begin
    Σ = [[0], [1], [-1]]
    @test Ada(Σ) == [[0], [1], [2]]
    @test Mat(Ada(Σ), Σ) == [1 1 1; 0 1 -1; 0 1 1]
    Σ = [[1], [-1]]
    @test Ada(Σ) == [[0], [1]]
    @test Mat(Ada(Σ), Σ) == [1 1; 1 -1]
    Σ = [[0], [1]]
    @test Ada(Σ) == [[0], [1]]
    @test Mat(Ada(Σ), Σ) == [1 1; 0 1]
    Σ = [[0], [-1]]
    @test Ada(Σ) == [[0], [1]]
    @test Mat(Ada(Σ), Σ) == [1 1; 0 -1]
    Σ = [[0]]
    @test Ada(Σ) == [[0]]
    @test Mat(Ada(Σ), Σ) == [1;;]
    Σ = [[1]]
    @test Ada(Σ) == [[0]]
    @test Mat(Ada(Σ), Σ) == [1;;]
    Σ = [[-1]]
    @test Ada(Σ) == [[0]]
    @test Mat(Ada(Σ), Σ) == [1;;]
end

# Example 10.81
@testset "Algorithms -> Adapted family 2" begin
    Σ = [
        [0, 1, 0, 0], [0, 1, 0, 1], [0, 1, 0, -1], [0, 1, 1, -1],
        [1, -1, 0, 0], [1, -1, 0, 1], [1, -1, 0, -1], [-1, 0, 0, 0],
        [-1, 0, 0, -1], [-1, 0, 1, 1], [-1, 0, 1, -1]
    ]
    @test Ada(Σ) == [
        [0, 0, 0, 0], [1, 0, 0, 0], [2, 0, 0, 0], [0, 0, 1, 0],
        [1, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 1], [2, 0, 0, 1],
        [0, 0, 1, 1], [0, 0, 0, 2], [1, 0, 0, 2]
    ]
end

@testset "Algorithms -> Adapted family 3" begin
    s = 10
    Σ = [Vector{Int}()]
    for i in 1:s
        Σ = Σ ∧ [1, -1]
    end
    A = Ada(Σ)
    M = matrix(QQ, Mat(A, Σ))
    @test size(M) == (2^s, 2^s)
    M = inv(M)
    @test M[end, end] == QQ(1 // 2^s)
end

@testset "Algorithms -> Sign determination 1" begin
    R, () = polynomial_ring(QQ, Symbol[])
    R′, (x1, x2) = polynomial_ring(fraction_field(R), [:x1, :x2])
    f = ParametricIdeal([x1^2 + x2^2 - 1, x1^2 - x2])
    g₁ = x1 + x2
    g₂ = x1 - x2
    S = AlgebraicSolving.sign_determination(f, [g₁, g₂])
    @test S == ([[1, 1], [-1, -1]], [1, 1])
end

@testset "Algorithms -> Sign determination 2" begin
    R, () = polynomial_ring(QQ, Symbol[])
    R′, (x1, x2, x3) = polynomial_ring(fraction_field(R), [:x1, :x2, :x3])
    f = ParametricIdeal([x1^2 - 1, x2^2 - 1, x3^3 - x3])
    S = AlgebraicSolving.sign_determination(f, [x1, x2, x3])
    @test length(S[1]) == 12
end
