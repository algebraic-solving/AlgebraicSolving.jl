
function Base.hash(diagram::Diagram, h::UInt)
    return hash(diagram.id, h)
end

function Base.:(==)(d1::Diagram, d2::Diagram)
    return d1.id == d2.id
end

function Base.show(io::IO, diagram::Diagram)
    print(diagram.edges)
end

# Function to print the diagram
function print_diagram(diagram::Diagram, space::Int=0)
    for edge in diagram.edges
        println(" " ^ space * " ", edge[1])
        print_diagram(edge[2], space + 2)
    end
end

# Define a function to create a new hash state
function new_hashstate()
    return HashState(Dict{Vector{Edge}, Diagram}(), 0)
end

function troncate!(children::Vector{Edge}, i::Int)
    if length(children) > i
        resize!(children,i)
    end
end

# Function to create a monomial divisibility diagram with hash consing
function make_diagram(children::Vector{Edge}, hashstate::HashState)
    i = 0
    last = empty_diagram
    for j in 1:length(children)
        current_edge = children[j][2]
        if current_edge != last
            last = current_edge
            children[i+1] = children[j]
            i += 1
        end
    end

    troncate!(children, i)
    key = get(hashstate.hashtable, children, empty_diagram)

    if key != empty_diagram
        return hashstate.hashtable[children]
    end

    diagram = Diagram(hashstate.counter, children)
    hashstate.hashtable[children] = diagram
    hashstate.counter += 1
    return diagram
end

# Function to insert a prefix in the diagram
function insert_aux(diagram::Diagram, m::Monomial, i::Int, hashstate::HashState)
    if isempty(m.exps) || i == 1
        return diagram
    end

    i -= 1
    a = m.exps[i]
    new_diagram = Edge[]
    inserted = false

    for edge in diagram.edges
        if edge[1] == a
            push!(new_diagram, (a,insert_aux(edge[2], m, i, hashstate)))
            inserted = true
        else
            push!(new_diagram, edge)
        end
    end
    if !inserted
        push!(new_diagram, (a,insert_aux(make_diagram(Edge[], hashstate),m,i,hashstate)))
    end
    return make_diagram(new_diagram, hashstate)
end

function insert(diagram::Diagram, m::Monomial, hashstate::HashState)
    return insert_aux(diagram, m, length(m.exps)+1, hashstate)
end


# Function to insert a monomial in the diagram with memoization
function insertion_aux(diagram::Diagram, m::Monomial, i::Int, memo::Dict{Tuple{Diagram,Int}, Diagram}, hashstate::HashState)
    if diagram == empty_diagram
        return insert_aux(make_diagram(Edge[], hashstate), m, i, hashstate)
    end
    if i == 1
         return diagram
    end
    if isempty(diagram.edges)
        return diagram
    end
    key = get(memo, (diagram, i), empty_diagram)
    if key != empty_diagram
        return key
    end

    i -= 1
    new_diagram = Edge[]
    is_sub_diagram = false
    last_sub_diagram = empty_diagram
    for edge in diagram.edges
        if edge[1] < m.exps[i]
            push!(new_diagram, edge)
            last_sub_diagram = edge[2]
        else
            if edge[1] == m.exps[i]
                is_sub_diagram = true
            end
            if edge[1] > m.exps[i] && !is_sub_diagram
                push!(new_diagram, (m.exps[i], insertion_aux(last_sub_diagram, m, i, memo, hashstate)))
                is_sub_diagram = true
            end

            push!(new_diagram, (edge[1], insertion_aux(edge[2], m, i, memo, hashstate)))
        end
    end

    if !is_sub_diagram
        push!(new_diagram, (m.exps[i], insertion_aux(last_sub_diagram, m, i, memo, hashstate)))
    end

    final_diagram = make_diagram(new_diagram, hashstate)
    memo[(diagram, i)] = final_diagram
    return final_diagram
end

function insertion(diagram::Diagram, m::Monomial, hashstate::HashState)
    return insertion_aux(diagram, m, length(m.exps)+1, Dict{Tuple{Diagram,Int}, Diagram}(), hashstate)
end


# Function to create the monomial divisibility diagram of a list of monomials
function create_diagram(monomial_list::Vector{<:Monomial}, hashstate::HashState)
    diagram = empty_diagram
    for m in monomial_list
        diagram = insertion(diagram, m, hashstate)
    end
    return diagram
end

# Function that determines the number of nodes in the monomial divisibility diagram
function size_of_diagram(diagram::Diagram)
    if isempty(diagram.edges) || diagram == empty_diagram
        return 0
    end
    return 1 + sum(size_of_diagram(sub_diagram[2]) for sub_diagram in diagram.edges)
end

# Function that determines the number of distinct nodes in the monomial divisibility diagram
function number_of_distinct_nodes(diagram::Diagram, mem::Dict{Diagram, Diagram} = Dict{Diagram,Diagram}())
    if isempty(diagram.edges) || diagram == empty_diagram
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

# Function that test if a monomial is represented by a monomial divisibility diagram
function is_in_diagram(m::Monomial{N}, diagram::Diagram) where N
    sub_diagram = diagram

    @inbounds for j in N:-1:1
        exp = m.exps[j]
        i = find_nearest_index(sub_diagram.edges, exp)
        if i == -1
            return false
        else
            sub_diagram = sub_diagram.edges[i][2]
        end
    end
    return true
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



@doc Markdown.doc"""
    hilbert_series_mdd(diagram::Diagram)

Compute the Hilbert series of a given monomial divisibility diagram.

**Notes**:
* This requires a monomial divisibility diagram of an ideal.

# Examples
```jldoctest
julia> using AlgebraicSolving

julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"]);

# This represents the monomials [x*y, x*z, y*z].
julia> list_of_monomials = [monomial(SVector{3}([0,1,1])), monomial(SVector{3}([1,0,1])), monomial(SVector{3}([1,1,0]))]

julia> diagram = create_diagram(list_of_monomials)

julia> hilbert_series_mdd(diagram)
(-2*t - 1)//(t - 1)
```
"""
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



@doc Markdown.doc"""
    hilbert_series_mdd(diagram::Diagram)

Compute the multivariate Hilbert series of a given monomial divisibility diagram.

**Notes**:
* This requires a monomial divisibility diagram of an ideal.

# Examples
```jldoctest
julia> using AlgebraicSolving

julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"]);

# This represents the monomials [x*y, x*z, y*z].
julia> list_of_monomials = [monomial(SVector{3}([0,1,1])), monomial(SVector{3}([1,0,1])), monomial(SVector{3}([1,1,0]))]

julia> diagram = create_diagram(list_of_monomials)

julia> multi_hilbert_series_mdd(diagram)
(-2*x1*x2*x3 + x1*x2 + x1*x3 + x2*x3 - 1)//(x1*x2*x3 - x1*x2 - x1*x3 + x1 - x2*x3 + x2 + x3 - 1)
```
"""
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