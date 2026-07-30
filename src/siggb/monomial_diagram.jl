# Function to hash a diagram.
function Base.hash(diagram::Diagram, h::UInt)
    return hash(objectid(diagram.edges), hash(diagram.id, h))
end

# Function to test equality between two diagrams in O(1).
function Base.:(==)(d1::Diagram, d2::Diagram)
    return d1.id == d2.id &&
           d1.index == d2.index &&
           d1.edges === d2.edges
end

function Base.show(io::IO, diagram::Diagram)
    if diagram == EMPTY_DIAGRAM
        print(io, "EmptyDiagram")
    elseif isempty(diagram.edges)
        print(io, "Terminal(", diagram.index, ")")
    else
        print(io, diagram.edges)
    end
end

# Function to print the diagram as an n-tree.
function print_diagram(diagram::Diagram, space::Int=0)
    if isempty(diagram.edges)
        println(" " ^ space * "Terminal(", diagram.index, ")")
        return
    end
    for edge in diagram.edges
        println(" " ^ space * " ", edge[1])
        print_diagram(edge[2], space + 2)
    end
end

# This function is to create a new hash state.
function new_hashstate()
    hashtable = Dict{Tuple{Vector{Edge}, Int}, Diagram}()
    leaf = Diagram(0, Edge[], 0)
    hashtable[(leaf.edges, leaf.index)] = leaf
    memo = Dict{Tuple{Diagram,Int}, Diagram}()
    return HashState(hashtable, leaf, memo, 1, 0, 0)
end

# Function to troncate a list of edges.
function truncate!(children::Vector{Edge}, i::Int)
    if length(children) > i
        resize!(children,i)
    end
end

@inline function push_strict_edge!(children::Vector{Edge}, edge::Edge)
    if isempty(children) || children[end][2] != edge[2]
        push!(children, edge)
    end
    return children
end

# Intern an already reduced list of edges.
function intern_diagram(children::Vector{Edge},
                        hashstate::HashState,
                        index::Int=0)
    @assert isempty(children) || iszero(index)
    key = (children, index)
    diagram = get(hashstate.hashtable, key, EMPTY_DIAGRAM)
    diagram != EMPTY_DIAGRAM && return diagram

    diagram = Diagram(hashstate.counter, children, index)
    hashstate.hashtable[key] = diagram
    hashstate.counter += 1
    return diagram
end

# Function to create a monomial divisibility diagram with hash consing
function make_diagram(children::Vector{Edge},
                      hashstate::HashState,
                      index::Int=0)
    i = 0

    # In order to keep strict inclusions we remove the edges that are the same keeping the leftmost one.
    last = EMPTY_DIAGRAM
    for j in 1:length(children)
        current_edge = children[j][2]
        if current_edge != last
            last = current_edge
            children[i+1] = children[j]
            i += 1
        end
    end

    truncate!(children, i)
    return intern_diagram(children, hashstate, index)
end

# Function to insert a prefix in the diagram
function insert_aux(diagram::Diagram,
                    m::Monomial,
                    i::Int,
                    hashstate::HashState,
                    index::Int=0)
    if isempty(m.exps) || i == 1
        return diagram
    end

    i -= 1
    a = m.exps[i]
    new_diagram = Edge[]
    sizehint!(new_diagram, length(diagram.edges) + 1)
    inserted = false

    for edge in diagram.edges
        if edge[1] == a
            push_strict_edge!(
                new_diagram,
                (a, insert_aux(edge[2], m, i, hashstate, index)),
            )
            inserted = true
        else
            push_strict_edge!(new_diagram, edge)
        end
    end
    if !inserted
        terminal = make_diagram(Edge[], hashstate, index)
        push_strict_edge!(
            new_diagram,
            (a, insert_aux(terminal, m, i, hashstate, index)),
        )
    end
    return intern_diagram(new_diagram, hashstate)
end

function insert(diagram::Diagram,
                m::Monomial,
                hashstate::HashState,
                index::Int=0)
    return insert_aux(diagram, m, length(m.exps)+1,
                      hashstate, index)
end


# Function to insert a monomial in the diagram with memoization
function insertion_aux(
    diagram::Diagram,
    m::Monomial,
    i::Int,
    memo::Dict{Tuple{Diagram,Int}, Diagram},
    hashstate::HashState,
)
    if diagram == EMPTY_DIAGRAM
        return insert_aux(hashstate.leaf, m, i, hashstate)
    end
    if i == 1
        return diagram
    end
    if isempty(diagram.edges)
        return diagram
    end
    level = i
    key = get(memo, (diagram, level), EMPTY_DIAGRAM)
    if key != EMPTY_DIAGRAM
        return key
    end

    i -= 1
    new_diagram = Edge[]
    sizehint!(new_diagram, length(diagram.edges) + 1)
    is_sub_diagram = false
    last_sub_diagram = EMPTY_DIAGRAM
    for edge in diagram.edges
        if edge[1] < m.exps[i]
            push_strict_edge!(new_diagram, edge)
            last_sub_diagram = edge[2]
        else
            if edge[1] == m.exps[i]
                is_sub_diagram = true
            end
            if edge[1] > m.exps[i] && !is_sub_diagram
                child = insertion_aux(
                    last_sub_diagram, m, i, memo, hashstate
                )
                push_strict_edge!(new_diagram, (m.exps[i], child))
                is_sub_diagram = true
            end

            child = insertion_aux(edge[2], m, i, memo, hashstate)
            push_strict_edge!(new_diagram, (edge[1], child))
        end
    end

    if !is_sub_diagram
        child = insertion_aux(
            last_sub_diagram, m, i, memo, hashstate
        )
        push_strict_edge!(new_diagram, (m.exps[i], child))
    end

    final_diagram = intern_diagram(new_diagram, hashstate)
    memo[(diagram, level)] = final_diagram
    return final_diagram
end

function insertion(diagram::Diagram, m::Monomial, hashstate::HashState)
    diagram != EMPTY_DIAGRAM && is_in_diagram(m, diagram) && return diagram
    empty!(hashstate.insertion_memo)
    return insertion_aux(
        diagram,
        m,
        length(m.exps) + 1,
        hashstate.insertion_memo,
        hashstate,
    )
end

@inline function best_terminal(diagram::Diagram,
                               candidate::Int,
                               is_better,
                               hashstate::HashState)
    if iszero(diagram.index) || is_better(candidate, diagram.index)
        return make_diagram(Edge[], hashstate, candidate)
    end
    return diagram
end

# Overlay the divisibility cone of `m` with a labelled terminal. Outside the
# cone the diagram is unchanged; inside it the supplied comparator chooses
# between the new label and the label already present.
function insertion_best_aux(
    diagram::Diagram,
    m::Monomial,
    i::Int,
    memo::Dict{Tuple{Diagram,Int}, Diagram},
    hashstate::HashState,
    candidate::Int,
    is_better,
)
    if diagram == EMPTY_DIAGRAM
        terminal = make_diagram(Edge[], hashstate, candidate)
        return insert_aux(terminal, m, i, hashstate, candidate)
    end
    if i == 1
        return best_terminal(diagram, candidate, is_better, hashstate)
    end

    level = i
    cached = get(memo, (diagram, level), EMPTY_DIAGRAM)
    cached != EMPTY_DIAGRAM && return cached

    i -= 1
    new_diagram = Edge[]
    sizehint!(new_diagram, length(diagram.edges) + 1)
    has_threshold = false
    previous = EMPTY_DIAGRAM

    for edge in diagram.edges
        if edge[1] < m.exps[i]
            push_strict_edge!(new_diagram, edge)
            previous = edge[2]
            continue
        end

        if edge[1] == m.exps[i]
            has_threshold = true
        elseif !has_threshold
            child = insertion_best_aux(
                previous, m, i, memo, hashstate, candidate, is_better
            )
            push_strict_edge!(new_diagram, (m.exps[i], child))
            has_threshold = true
        end

        child = insertion_best_aux(
            edge[2], m, i, memo, hashstate, candidate, is_better
        )
        push_strict_edge!(new_diagram, (edge[1], child))
    end

    if !has_threshold
        child = insertion_best_aux(
            previous, m, i, memo, hashstate, candidate, is_better
        )
        push_strict_edge!(new_diagram, (m.exps[i], child))
    end

    result = intern_diagram(new_diagram, hashstate)
    memo[(diagram, level)] = result
    return result
end

function insertion_best(diagram::Diagram,
                        m::Monomial,
                        hashstate::HashState,
                        candidate::Int,
                        is_better)
    candidate > 0 ||
        throw(ArgumentError("a labelled terminal must have a positive index"))
    empty!(hashstate.insertion_memo)
    return insertion_best_aux(
        diagram,
        m,
        length(m.exps) + 1,
        hashstate.insertion_memo,
        hashstate,
        candidate,
        is_better,
    )
end


# Function to create the monomial divisibility diagram of a list of monomials
function create_diagram(monomial_list::Vector{<:Monomial}, hashstate::HashState = new_hashstate())
    diagram = EMPTY_DIAGRAM
    for m in monomial_list
        diagram = insertion(diagram, m, hashstate)
    end
    return diagram
end

# Function that determines the number of nodes in the monomial divisibility diagram
function size_of_diagram(diagram::Diagram)
    if isempty(diagram.edges) || diagram == EMPTY_DIAGRAM
        return 0
    end
    return 1 + sum(size_of_diagram(sub_diagram[2]) for sub_diagram in diagram.edges)
end

# Function that determines the number of distinct nodes in the monomial divisibility diagram
function number_of_distinct_nodes(diagram::Diagram, mem::Dict{Diagram, Diagram} = Dict{Diagram,Diagram}())
    if isempty(diagram.edges) || diagram == EMPTY_DIAGRAM
        return 0
    end
    if haskey(mem, diagram)
        return 0
    end
    number = 0
    for sub_diagram in diagram.edges
        if !haskey(mem, sub_diagram[2])
            number += 1 + number_of_distinct_nodes(sub_diagram[2], mem)
        end
    end
    mem[diagram] = diagram
    return number
end

# Count all nodes reachable from several roots, including terminal nodes.
# EMPTY_DIAGRAM is only a sentinel and is therefore not counted.
function number_of_distinct_dag_nodes(roots::AbstractVector{Diagram})
    seen = Set{Diagram}()
    stack = Diagram[]
    for root in roots
        root != EMPTY_DIAGRAM && push!(stack, root)
    end

    while !isempty(stack)
        diagram = pop!(stack)
        diagram in seen && continue
        push!(seen, diagram)
        for (_, child) in diagram.edges
            child != EMPTY_DIAGRAM && push!(stack, child)
        end
    end
    return length(seen)
end

# Function that determines the largest s such that s <= exp, returns -1 otherwise
function find_nearest_index(sub_diagram::Vector{Tuple{Exp, Diagram}}, exp::Exp)
    j = -1
    for k in 1:length(sub_diagram)
        if sub_diagram[k][1] > exp
            return j
        end
        j = k
    end
    return j
end

# Follow the unique path associated with `m`. The Boolean diagrams end in
# Terminal(0), whereas labelled diagrams end in Terminal(id).
function lookup_diagram(m::Monomial{N}, diagram::Diagram) where N
    diagram == EMPTY_DIAGRAM && return false, 0
    sub_diagram = diagram

    @inbounds for j in N:-1:1
        exp = m.exps[j]
        i = find_nearest_index(sub_diagram.edges, exp)
        if i == -1
            return false, 0
        else
            sub_diagram = sub_diagram.edges[i][2]
        end
    end
    return true, sub_diagram.index
end

# Function that tests if a monomial is represented by a monomial diagram.
function is_in_diagram(m::Monomial, diagram::Diagram)
    found, _ = lookup_diagram(m, diagram)
    return found
end

function terminal_index(m::Monomial, diagram::Diagram)
    found, index = lookup_diagram(m, diagram)
    return found ? index : 0
end

# Function that converts a monomial::MPolyRingElem into an object of type Monomial{N}
function convert_to_monomial(m::MPolyRingElem, R::MPolyRing, ::Val{N}) where N
    exts = first(collect(exponent_vectors(R(m))))
    sv = SVector{N, Int}(exts)
    return monomial(sv)
end

# Function that tests if a divides b
function divides(a::Monomial, b::Monomial)
    return all(ai <= bi for (ai, bi) in zip(a.exps, b.exps))
end

# Function that generates a random list of monomials
function generate_random_vectors(r::Int, d::Int, n::Int)
    N = 0
    list_of_monomials = Monomial[]
    while N < n
        rd = monomial(SVector{r}(rand(Exp(0):Exp(d-1), r)))
        if all(!divides(mon, rd) for mon in list_of_monomials)
            push!(list_of_monomials, rd)
        end
        N += 1
    end

    return list_of_monomials
end

# Function that naively test if a monomial is in a monomial ideal
function naive_is_in_ideal(mon::Monomial, list::Vector{Monomial})
    for i in 1:length(list)
        if divides(list[i], mon)
            return true
        end
    end
    return false
end

# Function that returns the list of all leaves of a monomial divisibility diagram
function get_leaves(diagram::Diagram)
    if isempty(diagram.edges)
        return Vector{Exp}[Exp[]]
    end

    leaves = Vector{Exp}[]
    for (label, sub_diagram) in diagram.edges
        for path in get_leaves(sub_diagram)
            push!(leaves, vcat([label],path))
        end
    end
    return leaves
end

function depth(diagram::Diagram)::Int
    if isempty(diagram.edges)
        return 0
    end

    return 1 + depth(diagram.edges[1][2])
end

function hilbert_series_mdd_aux(R, i::Int, diagram::Diagram)
    t = gens(R)[1]

    if i == 1
        somme = R(0)
        for j in 0:(diagram.edges[1][1]-1)
            somme += t^j
        end
        return (R(1)-t) * somme
    end

    hilbert_series = R(1) - t^(diagram.edges[1][1])
    for j in 1:length(diagram.edges)
        if j == length(diagram.edges)
            hilbert_series += t^(diagram.edges[j][1]) * hilbert_series_mdd_aux(R, i-1, diagram.edges[j][2])
        else
            hilbert_series += (t^(diagram.edges[j][1]) - t^(diagram.edges[j+1][1])) * hilbert_series_mdd_aux(R, i-1, diagram.edges[j][2])
        end
    end

    return hilbert_series
end

function hilbert_series_mdd(diagram::Diagram)
    n = depth(diagram)
    R, t = polynomial_ring(ZZ, :t)
    K = fraction_field(R)
    return K(hilbert_series_mdd_aux(R, n, diagram)) / (K(1)-K(t))^n
end

function multi_hilbert_series_mdd_aux(R, i::Int, diagram::Diagram)
    vars = gens(R)
    if i == 1
        somme = R(0)
        for j in 0:(diagram.edges[1][1]-1)
            somme += vars[1]^j
        end
        return (R(1)-vars[1]) * somme
    end

    hilbert_series = 1 - vars[i]^(diagram.edges[1][1])
    for j in 1:length(diagram.edges)
        if j == length(diagram.edges)
            hilbert_series += vars[i]^(diagram.edges[j][1]) * multi_hilbert_series_mdd_aux(R, i-1, diagram.edges[j][2])
        else
            hilbert_series += (vars[i]^(diagram.edges[j][1]) - vars[i]^(diagram.edges[j+1][1])) * multi_hilbert_series_mdd_aux(R, i-1, diagram.edges[j][2])
        end
    end

    return hilbert_series
end

function multi_hilbert_series_mdd(diagram::Diagram)
    n = depth(diagram)
    R, vars = polynomial_ring(ZZ, ["x$i" for i in 1:n])
    K = fraction_field(R)
    prodn = K(1)
    for x in vars
        prodn *= K(R(1)-x)
    end
    return K(multi_hilbert_series_mdd_aux(R, n, diagram)) / prodn
end
