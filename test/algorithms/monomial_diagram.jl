using Nemo
using Random
using StaticArrays

@testset "Algorithms -> Monomial Diagram" begin
    R, (x1,x2,x3,x4,x5,x6) = polynomial_ring(GF(101), ["x1", "x2", "x3", "x4", "x5", "x6"], internal_ordering=:degrevlex)
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

    redundant = AlgebraicSolving.monomial(
        AlgebraicSolving.SVector{3, AlgebraicSolving.Exp}(1, 0, 0)
    )
    multiple = AlgebraicSolving.monomial(
        AlgebraicSolving.SVector{3, AlgebraicSolving.Exp}(2, 1, 0)
    )
    redundant_hashstate = AlgebraicSolving.new_hashstate()
    redundant_diagram = AlgebraicSolving.create_diagram([redundant], redundant_hashstate)
    counter_before = redundant_hashstate.counter
    memo = redundant_hashstate.insertion_memo
    @test AlgebraicSolving.insertion(redundant_diagram, multiple, redundant_hashstate) === redundant_diagram
    @test redundant_hashstate.counter == counter_before
    @test redundant_hashstate.insertion_memo === memo
    @test AlgebraicSolving.make_diagram(AlgebraicSolving.Edge[], redundant_hashstate) === redundant_hashstate.leaf

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

const CanonicalTestN = 3

canonical_test_monomial(exponents) = AlgebraicSolving.monomial(
    SVector{CanonicalTestN, AlgebraicSolving.Exp}(
        AlgebraicSolving.Exp.(exponents)
    )
)

function new_canonical_test_fixture(; nindices::Int=4)
    basis = AlgebraicSolving.new_basis(
        AlgebraicSolving.init_basis_size,
        AlgebraicSolving.init_syz_size,
        nindices,
        Val(CanonicalTestN);
        use_canonical_mdd=true,
    )
    basis.input_load = 0
    basis.input_size = nindices
    basis.basis_load = 0
    basis.basis_offset = 1

    ht = AlgebraicSolving.initialize_basis_hash_table(Val(CanonicalTestN))
    AlgebraicSolving.insert_in_hash_table!(
        ht, canonical_test_monomial((0, 0, 0))
    )
    AlgebraicSolving.insert_in_hash_table!(
        ht, canonical_test_monomial((12, 12, 12))
    )
    AlgebraicSolving.fill_divmask!(ht)
    return basis, ht
end

function initialize_canonical_test_slot!(
    basis,
    ht,
    slot::Int,
    idx::Int,
    signature_exponents,
    ratio_exponents,
)
    signature = (
        AlgebraicSolving.SigIndex(idx),
        canonical_test_monomial(signature_exponents),
    )
    basis.sigs[slot] = signature
    basis.sigmasks[slot] = (
        AlgebraicSolving.SigIndex(idx),
        AlgebraicSolving.divmask(
            AlgebraicSolving.monomial(signature), ht.divmap, ht.ndivbits
        ),
    )
    basis.sigratios[slot] = canonical_test_monomial(ratio_exponents)
    basis.lm_masks[slot] = zero(AlgebraicSolving.DivMask)
    basis.monomials[slot] = AlgebraicSolving.MonIdx[]
    basis.coefficients[slot] = AlgebraicSolving.Coeff[]
    basis.is_red[slot] = false
    basis.mod_rep_known[slot] = Bool[]
    basis.mod_reps[slot] = AlgebraicSolving.Polynomial[]
    basis.rewrite_nodes[slot + 1] = [-1, 1]
    return signature
end

function add_canonical_test_rule!(
    basis,
    ht,
    idx::Int,
    signature_exponents,
    ratio_exponents,
)
    slot = basis.basis_load + 1
    signature = initialize_canonical_test_slot!(
        basis,
        ht,
        slot,
        idx,
        signature_exponents,
        ratio_exponents,
    )
    basis.basis_load = slot
    id = AlgebraicSolving.register_basis_id!(basis, slot)
    AlgebraicSolving.update_canonical_diagram!(basis, signature, id)
    return slot, id
end

function canonical_test_winners(basis, ht, idx::Int, exponents)
    signature = (
        AlgebraicSolving.SigIndex(idx),
        canonical_test_monomial(exponents),
    )
    signature_mask = AlgebraicSolving.divmask(
        AlgebraicSolving.monomial(signature), ht.divmap, ht.ndivbits
    )
    return (
        scan=AlgebraicSolving.find_canonical_rewriter_rat_scan(
            basis, signature, signature_mask
        ),
        mdd=AlgebraicSolving.find_canonical_rewriter_rat_mdd(
            basis, signature, signature_mask
        ),
    )
end

@testset "Algorithms -> Indexed monomial diagram" begin
    bool_state = AlgebraicSolving.new_hashstate()
    x = canonical_test_monomial((1, 0, 0))
    xy = canonical_test_monomial((1, 1, 0))
    bool_diagram = AlgebraicSolving.create_diagram([x, xy], bool_state)

    found, terminal = AlgebraicSolving.lookup_diagram(xy, bool_diagram)
    @test found
    @test terminal == 0
    @test bool_state.leaf.index == 0
    @test bool_state.leaf != AlgebraicSolving.EMPTY_DIAGRAM
    @test bool_state.leaf != AlgebraicSolving.new_hashstate().leaf
    @test AlgebraicSolving.make_diagram(
        AlgebraicSolving.Edge[], bool_state
    ) === bool_state.leaf

    # lm_diagram and koszul_diagram are ordinary Boolean diagrams and share
    # Terminal(0) when they use the same hash-consing state.
    basis = AlgebraicSolving.new_basis(
        AlgebraicSolving.init_basis_size,
        AlgebraicSolving.init_syz_size,
        1,
        Val(CanonicalTestN),
    )
    basis.lm_diagram = AlgebraicSolving.insertion(
        basis.lm_diagram, x, basis.hashstate
    )
    basis.koszul_diagram = AlgebraicSolving.insertion(
        basis.koszul_diagram, x, basis.hashstate
    )
    @test basis.lm_diagram === basis.koszul_diagram
    @test AlgebraicSolving.lookup_diagram(
        x, basis.lm_diagram
    ) == (true, 0)

    # A = xy, B = y². A has the better priority and must therefore win on
    # their overlap xy².
    priority = [1, 2]
    better = (candidate, current) ->
        priority[candidate] < priority[current] ||
        (
            priority[candidate] == priority[current] &&
            candidate > current
        )
    labelled_state = AlgebraicSolving.new_hashstate()
    labelled = AlgebraicSolving.insertion_best(
        AlgebraicSolving.EMPTY_DIAGRAM,
        canonical_test_monomial((1, 1, 0)),
        labelled_state,
        1,
        better,
    )
    labelled = AlgebraicSolving.insertion_best(
        labelled,
        canonical_test_monomial((0, 2, 0)),
        labelled_state,
        2,
        better,
    )
    @test AlgebraicSolving.terminal_index(
        canonical_test_monomial((1, 2, 0)), labelled
    ) == 1
    @test AlgebraicSolving.terminal_index(
        canonical_test_monomial((0, 2, 0)), labelled
    ) == 2
    @test AlgebraicSolving.terminal_index(
        canonical_test_monomial((0, 1, 0)), labelled
    ) == 0
    labelled_size =
        AlgebraicSolving.number_of_distinct_dag_nodes([labelled])
    @test labelled_size > 0
    @test AlgebraicSolving.number_of_distinct_dag_nodes(
        [labelled, labelled, AlgebraicSolving.EMPTY_DIAGRAM]
    ) == labelled_size
    @test AlgebraicSolving.number_of_distinct_dag_nodes(
        [labelled_state.leaf]
    ) == 1

    # Equal priorities reproduce the historical "last inserted wins" rule.
    tie_priority = [1, 1]
    tie_better = (candidate, current) ->
        tie_priority[candidate] < tie_priority[current] ||
        (
            tie_priority[candidate] == tie_priority[current] &&
            candidate > current
        )
    tie_state = AlgebraicSolving.new_hashstate()
    tie_diagram = AlgebraicSolving.insertion_best(
        AlgebraicSolving.EMPTY_DIAGRAM,
        canonical_test_monomial((1, 0, 0)),
        tie_state,
        1,
        tie_better,
    )
    tie_diagram = AlgebraicSolving.insertion_best(
        tie_diagram,
        canonical_test_monomial((0, 1, 0)),
        tie_state,
        2,
        tie_better,
    )
    @test AlgebraicSolving.terminal_index(
        canonical_test_monomial((1, 1, 0)), tie_diagram
    ) == 2
end

@testset "Algorithms -> Canonical rewriter MDD" begin
    basis, ht = new_canonical_test_fixture()
    slot_A, _ = add_canonical_test_rule!(
        basis, ht, 1, (1, 1, 0), (0, 0, 0)
    )
    slot_B, _ = add_canonical_test_rule!(
        basis, ht, 1, (0, 2, 0), (1, 0, 0)
    )
    for (query, expected) in (
        ((1, 1, 0), slot_A),
        ((0, 2, 0), slot_B),
        ((1, 2, 0), slot_A),
        ((4, 5, 0), slot_A),
        ((0, 1, 0), 0),
    )
        @test canonical_test_winners(basis, ht, 1, query) ==
              (scan=expected, mdd=expected)
    end

    # The signature index is a categorical first level.
    @test canonical_test_winners(basis, ht, 2, (1, 2, 0)) ==
          (scan=0, mdd=0)
    @test canonical_test_winners(basis, ht, 0, (1, 2, 0)) ==
          (scan=0, mdd=0)
    slot_index_two, _ = add_canonical_test_rule!(
        basis, ht, 2, (0, 0, 0), (0, 0, 0)
    )
    @test canonical_test_winners(basis, ht, 2, (1, 2, 0)) ==
          (scan=slot_index_two, mdd=slot_index_two)
    @test canonical_test_winners(basis, ht, 1, (1, 2, 0)) ==
          (scan=slot_A, mdd=slot_A)

    # Equal sigratios: the later stable ID wins on the overlap.
    tie_basis, tie_ht = new_canonical_test_fixture()
    slot_x, _ = add_canonical_test_rule!(
        tie_basis, tie_ht, 1, (1, 0, 0), (0, 0, 0)
    )
    slot_y, _ = add_canonical_test_rule!(
        tie_basis, tie_ht, 1, (0, 1, 0), (0, 0, 0)
    )
    @test canonical_test_winners(tie_basis, tie_ht, 1, (1, 0, 0)) ==
          (scan=slot_x, mdd=slot_x)
    @test canonical_test_winners(tie_basis, tie_ht, 1, (1, 1, 0)) ==
          (scan=slot_y, mdd=slot_y)

    # Three overlapping cones exercise all changes of canonical region.
    three_basis, three_ht = new_canonical_test_fixture()
    slot_x2, _ = add_canonical_test_rule!(
        three_basis, three_ht, 1, (2, 0, 0), (2, 0, 0)
    )
    slot_y2, _ = add_canonical_test_rule!(
        three_basis, three_ht, 1, (0, 2, 0), (0, 0, 0)
    )
    slot_xy, _ = add_canonical_test_rule!(
        three_basis, three_ht, 1, (1, 1, 0), (1, 0, 0)
    )
    for (query, expected) in (
        ((2, 0, 0), slot_x2),
        ((0, 2, 0), slot_y2),
        ((1, 1, 0), slot_xy),
        ((2, 1, 0), slot_xy),
        ((1, 2, 0), slot_y2),
        ((2, 2, 0), slot_y2),
    )
        @test canonical_test_winners(three_basis, three_ht, 1, query) ==
              (scan=expected, mdd=expected)
    end

    # Differential test against the untouched historical scan.
    random_basis, random_ht = new_canonical_test_fixture()
    rng = MersenneTwister(0x6d6464)
    for _ in 1:200
        add_canonical_test_rule!(
            random_basis,
            random_ht,
            rand(rng, 1:4),
            ntuple(_ -> rand(rng, 0:6), CanonicalTestN),
            ntuple(_ -> rand(rng, -2:5), CanonicalTestN),
        )
    end
    for _ in 1:300
        idx = rand(rng, 1:4)
        query = ntuple(_ -> rand(rng, 0:10), CanonicalTestN)
        winners = canonical_test_winners(
            random_basis, random_ht, idx, query
        )
        @test winners.mdd == winners.scan
    end
end

@testset "Algorithms -> Stable canonical IDs" begin
    # Rebuilding after is_red reveals the previous second-best rewriter.
    basis, ht = new_canonical_test_fixture()
    first_slot, _ = add_canonical_test_rule!(
        basis, ht, 1, (1, 0, 0), (0, 0, 0)
    )
    second_slot, _ = add_canonical_test_rule!(
        basis, ht, 1, (0, 1, 0), (0, 0, 0)
    )
    @test canonical_test_winners(basis, ht, 1, (1, 1, 0)) ==
          (scan=second_slot, mdd=second_slot)
    old_state = basis.canonical_hashstate
    basis.is_red[second_slot] = true
    AlgebraicSolving.rebuild_canonical_diagrams!(basis)
    @test basis.canonical_hashstate !== old_state
    @test canonical_test_winners(basis, ht, 1, (1, 1, 0)) ==
          (scan=first_slot, mdd=first_slot)

    # Garbage collection invalidates deleted IDs, updates moved IDs, rebuilds
    # the roots, and never recycles an ID.
    gc_basis, gc_ht = new_canonical_test_fixture()
    _, first_id = add_canonical_test_rule!(
        gc_basis, gc_ht, 1, (1, 0, 0), (2, 0, 0)
    )
    deleted_slot, deleted_id = add_canonical_test_rule!(
        gc_basis, gc_ht, 1, (0, 1, 0), (1, 0, 0)
    )
    moved_slot, moved_id = add_canonical_test_rule!(
        gc_basis, gc_ht, 1, (0, 0, 1), (0, 0, 0)
    )
    pairset = AlgebraicSolving.init_pairset(Val(CanonicalTestN))
    AlgebraicSolving.garbage_collect!(
        gc_basis,
        pairset,
        AlgebraicSolving.NoTracer(),
        [deleted_slot],
    )
    @test AlgebraicSolving.slot_from_id(gc_basis, deleted_id) == 0
    @test AlgebraicSolving.slot_from_id(gc_basis, moved_id) ==
          moved_slot - 1
    @test canonical_test_winners(gc_basis, gc_ht, 1, (1, 1, 1)).mdd ==
          canonical_test_winners(gc_basis, gc_ht, 1, (1, 1, 1)).scan
    _, new_id = add_canonical_test_rule!(
        gc_basis, gc_ht, 1, (2, 2, 2), (3, 0, 0)
    )
    @test new_id > max(first_id, deleted_id, moved_id)

    # A slot shift leaves terminal labels untouched and changes only id→slot.
    shift_basis = AlgebraicSolving.new_basis(
        AlgebraicSolving.init_basis_size,
        AlgebraicSolving.init_syz_size,
        1,
        Val(CanonicalTestN);
        use_canonical_mdd=true,
    )
    shift_ht = gc_ht
    shift_basis.input_load = 1
    shift_basis.input_size = 1
    shift_basis.basis_load = 1
    shift_basis.basis_offset = 2
    input_signature = initialize_canonical_test_slot!(
        shift_basis, shift_ht, 1, 1, (0, 0, 0), (0, 0, 0)
    )
    input_id = AlgebraicSolving.register_basis_id!(shift_basis, 1)
    AlgebraicSolving.update_canonical_diagram!(
        shift_basis, input_signature, input_id
    )
    computed_slot, computed_id = add_canonical_test_rule!(
        shift_basis, shift_ht, 1, (1, 0, 0), (0, 0, 0)
    )
    tracer = AlgebraicSolving.new_tracer()
    tracer.basis_ind_to_mat[computed_slot] = 0
    empty_pairset = AlgebraicSolving.init_pairset(Val(CanonicalTestN))
    old_root = shift_basis.canonical_diagrams[1]
    AlgebraicSolving.make_room_new_input_el!(
        shift_basis, empty_pairset, tracer
    )
    shifted_slot = computed_slot + shift_basis.input_size
    @test AlgebraicSolving.slot_from_id(shift_basis, computed_id) ==
          shifted_slot
    @test shift_basis.canonical_diagrams[1] === old_root
    @test canonical_test_winners(
        shift_basis, shift_ht, 1, (1, 0, 0)
    ) == (scan=shifted_slot, mdd=shifted_slot)
end

@testset "Algorithms -> Canonical scan/MDD end-to-end" begin
    R, (x, y, z) = polynomial_ring(
        GF(17), ["x", "y", "z"], internal_ordering=:degrevlex
    )
    systems = Vector{Vector{typeof(x)}}()
    push!(systems, [x^2 + y*z, x*y + z^2, y^2 + x*z])

    Rk, _ = polynomial_ring(
        GF(17), ["x$i" for i in 0:3], internal_ordering=:degrevlex
    )
    push!(systems, katsura(Rk).gens)

    Rc, _ = polynomial_ring(
        GF(17), ["x$i" for i in 1:4], internal_ordering=:degrevlex
    )
    push!(systems, cyclic(Rc).gens)

    for system in systems
        scan = AlgebraicSolving._sig_groebner_basis(
            system; use_canonical_mdd=false
        )
        mdd = AlgebraicSolving._sig_groebner_basis(
            system; use_canonical_mdd=true
        )
        @test mdd == scan
    end
end
