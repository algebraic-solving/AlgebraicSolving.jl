using Nemo
using Random
using StaticArrays

@testset "Algorithms -> Staircase trees" begin
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
    staircase_tree1 = AlgebraicSolving.create_staircase_tree(list_of_monomials1, hashstate)
    @test AlgebraicSolving.size_of_tree(staircase_tree1) == 155 
    @test AlgebraicSolving.number_of_distinct_trees(staircase_tree1) == 21

    hashstate = AlgebraicSolving.new_hashstate()
    lst = AlgebraicSolving.generate_random_vectors(10, 5, 10)
    
    T1 = AlgebraicSolving.create_staircase_tree(lst, hashstate)
    for k in 1:100
        lst2 = shuffle(lst)
        T2 = AlgebraicSolving.create_staircase_tree(lst2, hashstate)
        @test T1 == T2
    end

    for k in 1:100
        mon = AlgebraicSolving.generate_random_vectors(10, 5, 1)[1]
        list_reversed = []
        for m in mon.exps
            push!(list_reversed, m) 
        end
        tree_mon = AlgebraicSolving.monomial(SVector{10}(reverse(list_reversed)))
        @test AlgebraicSolving.naive_is_in_ideal(mon, lst) == AlgebraicSolving.is_in_tree(tree_mon, T1)
    end

    leaves = AlgebraicSolving.get_leaves(T1)
    @test is_empty([leaf for leaf in leaves if !AlgebraicSolving.naive_is_in_ideal(AlgebraicSolving.monomial(AlgebraicSolving.SVector{10}(reverse(leaf))), lst)])
end