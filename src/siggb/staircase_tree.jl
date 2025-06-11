# Change monomial encoding: Vector{Int} to Monomial{N}

const Exp = Int16

# Define a Tree structure to represent the staircase tree
struct Tree
    id::Int
    edges::Vector{Tuple{Exp, Tree}}
end

function Base.:(==)(t1::Tree, t2::Tree)
    return t1.id == t2.id
end

function Base.show(io::IO, tree::Tree)
    print(tree.edges)
end

# Define a constant for an empty tree
const nothing_tree = Tree(-1, Tuple{Int,Tree}[])

# Function to print the tree structure
function print_tree(tree::Tree, space::Int=0)
    for edge in tree.edges
        println(" " ^ space * " ", edge[1])
        print_tree(edge[2], space + 2)
    end
end

# Define an Edge as a tuple of an exponent and a Tree
const Edge = Tuple{Exp, Tree}

# Define a mutable struct to hold the hash state
mutable struct HashState
    hashtable::Dict{Vector{Edge}, Tree}
    counter::Int
end

# Define a function to create a new hash state
function new_hashstate()
    return HashState(Dict{Vector{Edge}, Tree}(), 0)
end

# Function to create a trees with hash consing
function make_Tree(children::Vector{Edge}, hashstate::HashState)
    edges = []
    last = nothing
    for child in children
        if isnothing(last) || child[2] != last
            push!(edges, child)
            last = child[2]
        end
    end

    if haskey(hashstate.hashtable, edges)
        return hashstate.hashtable[edges]
    end

    tree = Tree(hashstate.counter, edges)
    hashstate.counter += 1
    hashstate.hashtable[edges] = tree
    return tree
end

# Function to insert a prefix in the tree
function insert_aux(tree::Tree, m::Vector{Exp}, i::Int, hashstate::HashState)
    if m == [] || i == 1
        return tree
    end

    i -= 1
    a = m[i]
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
        push!(newtree, (a,insert_aux(make_Tree(Edge[], hashstate),m,i,hashstate)))
    end
    return make_Tree(newtree, hashstate)
end

function insert(tree::Tree, m::Vector{Exp}, hashstate::HashState)
    return insert_aux(tree, m, length(m)+1, hashstate)
end


# Function to insert a monomial into the tree with memoization
function insertion_aux(tree::Tree, m::Vector{Exp}, i::Int, memo::Dict{Tuple{Tree,Int}, Tree}, hashstate::HashState)
    if tree == nothing_tree
        return insert_aux(make_Tree(Edge[], hashstate), m, i, hashstate)
    end
    if i == 1
        return tree
    end
    if tree == make_Tree(Edge[], hashstate)
        return tree
    end
    if haskey(memo, (tree, i))
        return memo[(tree, i)]
    end

    i -= 1
    newtree = Edge[]
    is_subtree = false
    last_subtree = nothing_tree
    for edge in tree.edges
        if edge[1] < m[i]
            push!(newtree, edge)
            last_subtree = edge[2]
        else
            if edge[1] == m[i]
                is_subtree = true
            end
            push!(newtree, (edge[1], insertion_aux(edge[2], m, i, memo, hashstate)))
        end
    end
            
    if !is_subtree
        push!(newtree, (m[i], insertion_aux(last_subtree, m, i, memo, hashstate)))
        newtree = sort(newtree, by = x -> x[1])
    end
    final_tree = make_Tree(newtree, hashstate)
    memo[(tree, i)] = final_tree
    return final_tree
end

function insertion(tree::Tree, m::Vector{Exp}, hashstate::HashState)
    return insertion_aux(tree, m, length(m)+1, Dict{Tuple{Tree,Int}, Tree}(), hashstate)
end


# Function to create the staircase tree of a list of monomials
function creation(monomial_list::Vector{Vector{Exp}})
    tree = nothing_tree
    hashstate = new_hashstate()
    for m in monomial_list
        tree = insertion(tree, m, hashstate)
    end
    return tree
end


T = creation([Int16[1,1,1], Int16[0,2,2]])
print_tree(T)