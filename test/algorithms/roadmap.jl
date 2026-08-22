
# Helpers

"""
    validate_roadmap_structure(node::RMnode, expected_depth::Int)

Recursively traverses the Roadmap tree to ensure its structural and
algebraic integrity according to the algorithm's definition.
"""
function validate_roadmap_structure(node, expected_depth=0)
    # The base point length must exactly match the depth in the tree
    @test length(node.base_pt) == expected_depth

    for child in node.children
        # Children must strictly extend the parent's base_pt
        @test length(child.base_pt) == expected_depth + 1
        # All children must share the same polar equations (only differ by fiber)
        @test child.polar_eqs == node.children[1].polar_eqs

        # Parent's exact sequence is preserved in the child
        if expected_depth > 0
            @test child.base_pt[1:expected_depth] == node.base_pt
        end

        validate_roadmap_structure(child, expected_depth + 1)
    end
end

"""
    max_tree_depth(node::RMnode)

Computes the maximum depth of the roadmap tree to ensure recursion is happening.
"""
function max_tree_depth(node)
    isempty(node.children) && return 0
    return 1 + maximum(max_tree_depth(c) for c in node.children)
end

# Tests

@testset "Algorithms -> Roadmap" begin

    R3, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])

    ### Edge Cases: Low Dimension Limits

    # Case 1: 0-dimensional variety (a single point)
    I_0 = AlgebraicSolving.Ideal([x, y, z])
    C_0 = AlgebraicSolving.Ideal([x, y, z]) # Query point is the origin

    rm_0 = roadmap(I_0, C_0)

    @test rm_0.initial_ideal == I_0
    @test isempty(rm_0.root.children)
    @test length(rm_0.root.base_pt) == 0

    # Case 2: 1-dimensional variety (a curve in 3D)
    I_1 = AlgebraicSolving.Ideal([x^2 + y^2 - 1, z]) # A circle in the xy-plane
    C_1 = AlgebraicSolving.Ideal([x - 1, y, z])      # Query point (1,0,0)

    rm_1 = roadmap(I_1, C_1)
    @test isempty(rm_1.root.children)
    @test rm_1.root.polar_eqs == I_1.gens

    ### Standard Case: 3D Sphere

    # V: x^2 + y^2 + z^2 - 1 = 0
    I_sphere = AlgebraicSolving.Ideal([x^2 + y^2 + z^2 - 1])

    # P: North Pole (0,0,1) and South Pole (0,0,-1)
    # Defined by the ideal <x, y, z^2 - 1>
    C_sphere = AlgebraicSolving.Ideal([x, y, z^2 - 1])

    rm_sphere = roadmap(I_sphere, C_sphere)

    # Validate root properties
    @test rm_sphere.initial_ideal == I_sphere
    validate_roadmap_structure(rm_sphere.root, 0)

    @test !isempty(rm_sphere.root.children)
    @test max_tree_depth(rm_sphere.root) == 1

    rm_eqs = all_eqs(rm_sphere)

    @test nb_nodes(rm_sphere) == length(rm_eqs)
    @test all( isone(dimension(node)) for node in rm_eqs)

    ## Limit Case: Empty Real Variety

    # V: x^2 + y^2 + z^2 + 1 = 0 (No real solutions)
    I_empty = AlgebraicSolving.Ideal([x^2 + y^2 + z^2 + 1])
    C_empty = AlgebraicSolving.Ideal([x, y, z^2 + 1])

    # The algorithm should run without throwing bounds errors
    try
        rm_empty = roadmap(I_empty, C_empty)
        @test rm_empty.initial_ideal == I_empty
        validate_roadmap_structure(rm_empty.root, 0)
    catch e
        @test false # Fail test if the empty variety throws an unexpected exception
    end

    # Higher Dimension: 4 Variables
    R4, (x1, x2, x3, x4) = polynomial_ring(QQ, ["x1", "x2", "x3", "x4"])

    # V: 3-Sphere in 4D space
    I_4 = AlgebraicSolving.Ideal([x1^2 + x2^2 + x3^2 + x4^2 - 2])
    C_4 = AlgebraicSolving.Ideal([x1 - 1, x2 - 2, x3, x4])

    rm_4 = roadmap(I_4, C_4)
    validate_roadmap_structure(rm_4.root, 0)

    # Because we start at dim 3, the recursion must go deeper than the 3D cases
    @test max_tree_depth(rm_4.root) >= 2

    # Empty Query Points
    I_sphere = AlgebraicSolving.Ideal([x^2 + y^2 + z^2 - 1])
    C_empty = AlgebraicSolving.Ideal([R3(1)])

    rm_no_queries = roadmap(I_sphere, C_empty)
    validate_roadmap_structure(rm_no_queries.root, 0)
    @test !isempty(rm_no_queries.root.children)

end