
function Base.hash(tree::Tree, h::UInt)
    return hash(tree.id, h)
end

function Base.:(==)(t1::Tree, t2::Tree)
    return t1.id == t2.id
end

function Base.show(io::IO, tree::Tree)
    print(tree.edges)
end

# Function to print the tree structure
function print_tree(tree::Tree, space::Int=0)
    for edge in tree.edges
        println(" " ^ space * " ", edge[1])
        print_tree(edge[2], space + 2)
    end
end

# Define a function to create a new hash state
function new_hashstate()
    return HashState(Dict{Vector{Edge}, Tree}(), 0)
end

function troncate!(children::Vector{Edge}, i::Int)
    if length(children) > i
        resize!(children,i)
    end
end

# Function to create a trees with hash consing
function make_tree(children::Vector{Edge}, hashstate::HashState)
    i = 0
    last = nothing_tree
    for j in 1:length(children)
        current_edge = children[j][2]
        if current_edge != last
            last = current_edge
            children[i+1] = children[j]
            i += 1
        end
    end

    troncate!(children, i)
    key = get(hashstate.hashtable, children, nothing_tree)

    if key != nothing_tree
        return hashstate.hashtable[children]
    end

    tree = Tree(hashstate.counter, children)
    hashstate.hashtable[children] = tree
    hashstate.counter += 1
    return tree
end

# Function to insert a prefix in the tree
function insert_aux(tree::Tree, m::Monomial, i::Int, hashstate::HashState)
    if isempty(m.exps) || i == 1
        return tree
    end

    i -= 1
    a = m.exps[i]
    newtree = Edge[]
    inserted = false

    for edge in tree.edges
        if edge[1] == a
            push!(newtree, (a,insert_aux(edge[2], m, i, hashstate)))
            inserted = true
        else
            push!(newtree, edge)
        end
    end
    if !inserted
        push!(newtree, (a,insert_aux(make_tree(Edge[], hashstate),m,i,hashstate)))
    end
    return make_tree(newtree, hashstate)
end

function insert(tree::Tree, m::Monomial, hashstate::HashState)
    return insert_aux(tree, m, length(m.exps)+1, hashstate)
end


# Function to insert a monomial into the tree with memoization
function insertion_aux(tree::Tree, m::Monomial, i::Int, memo::Dict{Tuple{Tree,Int}, Tree}, hashstate::HashState)
    if tree == nothing_tree
        return insert_aux(make_tree(Edge[], hashstate), m, i, hashstate)
    end
    if i == 1
         return tree
    end
    if isempty(tree.edges)
        return tree
    end
    key = get(memo, (tree, i), nothing_tree)
    if key != nothing_tree
        return key
    end

    i -= 1
    newtree = Edge[]
    is_subtree = false
    last_subtree = nothing_tree
    for edge in tree.edges
        if edge[1] < m.exps[i]
            push!(newtree, edge)
            last_subtree = edge[2]
        else
            if edge[1] == m.exps[i]
                is_subtree = true
            end
            if edge[1] > m.exps[i] && !is_subtree
                push!(newtree, (m.exps[i], insertion_aux(last_subtree, m, i, memo, hashstate)))
                is_subtree = true
            end

            push!(newtree, (edge[1], insertion_aux(edge[2], m, i, memo, hashstate)))
        end
    end

    if !is_subtree
        push!(newtree, (m.exps[i], insertion_aux(last_subtree, m, i, memo, hashstate)))
    end

    final_tree = make_tree(newtree, hashstate)
    memo[(tree, i)] = final_tree
    return final_tree
end

function insertion(tree::Tree, m::Monomial, hashstate::HashState)
    return insertion_aux(tree, m, length(m.exps)+1, Dict{Tuple{Tree,Int}, Tree}(), hashstate)
end


# Function to create the staircase tree of a list of monomials
function create_staircase_tree(monomial_list::Vector{<:Monomial}, hashstate::HashState)
    tree = nothing_tree
    for m in monomial_list
        tree = insertion(tree, m, hashstate)
    end
    return tree
end

# Function that determines the number of nodes in the staircase tree
function size_of_tree(tree::Tree)
    if isempty(tree.edges) || tree == nothing_tree
        return 0
    end
    return 1 + sum(size_of_tree(subtree[2]) for subtree in tree.edges)
end

# Function that determines the number of distinct subtrees in the staircase tree
function number_of_distinct_trees(tree::Tree, mem::Dict{Tree, Tree} = Dict{Tree,Tree}())
    if isempty(tree.edges) || tree == nothing_tree
        return 0
    end
    if haskey(mem, tree)
        return 0
    end
    number = 0
    for subtree in tree.edges
        if !haskey(mem, subtree[2])
            number += 1 + number_of_distinct_trees(subtree[2], mem)
        end
    end
    mem[tree] = tree
    return number
end

# Function that determines the largest s such that s <= exp, returns -1 otherwise
function find_nearest_index(subtree::Vector{Tuple{Exp, Tree}}, exp::Exp)
    j = -1
    for k in 1:length(subtree)
        if subtree[k][1] > exp
            return j
        end
        j = k
    end
    return j
end

# Function that test if a monomial is represented by a staircase tree
function is_in_tree(m::Monomial{N}, tree::Tree) where N
    subtree = tree

    @inbounds for j in N:-1:1
        exp = m.exps[j]
        i = find_nearest_index(subtree.edges, exp)
        if i == -1
            return false
        else
            subtree = subtree.edges[i][2]
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

# Function that returns the list of all leaves of tree
function get_leaves(tree::Tree)
    if isempty(tree.edges)
        return Vector{Exp}[Exp[]]
    end

    leaves = Vector{Exp}[]
    for (label, subtree) in tree.edges
        for path in get_leaves(subtree)
            push!(leaves, vcat([label],path))
        end
    end
    return leaves
end
