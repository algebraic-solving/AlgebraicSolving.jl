# Change monomial encoding: Vector{Int} to Monomial{N}

struct Tree
    id::Int
    edges::Vector{Tuple{Exp, Tree}}
end

const Edge = Tuple{Int, Tree}

mutable struct HashState
    hashtable::Dict{Vector{Edge}, Tree}
    counter::Int
end

function new_hashstate()
    return HashState(Dict{Vector{Edge}, Tree}(), 0)
end
