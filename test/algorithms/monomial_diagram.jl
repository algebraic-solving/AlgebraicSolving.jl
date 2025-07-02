using Nemo
using Random
using StaticArrays

@testset "Algorithms -> Monomial Diagram" begin
    R, (x1,x2,x3,x4,x5,x6) = polynomial_ring(QQ, ["x1", "x2", "x3", "x4", "x5", "x6"], internal_ordering=:degrevlex)
    list_of_polynomials1 = [x1+2*x2+2*x3+2*x4+2*x5+2*x6-1, 
    x1^2+2*x2^2+2*x3^2+2*x4^2+2*x5^2+2*x6^2-x1,
    2*x1*x2+2*x2*x3+2*x3*x4+2*x4*x5+2*x5*x6-x2,
    x2^2+2*x1*x3+2*x2*x4+2*x3*x5+2*x4*x6-x3,
    2*x2*x3+2*x1*x4+2*x2*x5+2*x3*x6-x4,
    x3^2+2*x2*x4+2*x1*x5+2*x2*x6-x5]
    I = AlgebraicSolving.Ideal(list_of_polynomials1)
    G = AlgebraicSolving.groebner_basis(I)
    list_of_monomials1 = [AlgebraicSolving.convert_to_monomial(Nemo.leading_monomial(g), R, Val(6)) for g in G]
    hashstate = AlgebraicSolving.new_hashstate()
    diagram1 = AlgebraicSolving.create_diagram(list_of_monomials1, hashstate)
    @test AlgebraicSolving.size_of_diagram(diagram1) == 155 
    @test AlgebraicSolving.number_of_distinct_nodes(diagram1) == 21

    hashstate = AlgebraicSolving.new_hashstate()
    lst = AlgebraicSolving.generate_random_vectors(10, 5, 10)
    
    d1 = AlgebraicSolving.create_diagram(lst, hashstate)
    for k in 1:30
        lst2 = shuffle(lst)
        d2 = AlgebraicSolving.create_diagram(lst2, hashstate)
        @test d1 == d2
    end

    for k in 1:30
        mon = AlgebraicSolving.generate_random_vectors(10, 5, 1)[1]
        @test AlgebraicSolving.naive_is_in_ideal(mon, lst) == AlgebraicSolving.is_in_diagram(mon, d1)
    end

    leaves = AlgebraicSolving.get_leaves(d1)
    @test is_empty([leaf for leaf in leaves if !AlgebraicSolving.naive_is_in_ideal(AlgebraicSolving.monomial(AlgebraicSolving.SVector{10}(reverse(leaf))), lst)])
end