# Represents the isolated points and valid connections on a single box edge
struct BoxEdge{T}
    points::Vector{Vector{T}}
    indices_inside::Vector{Int}
end

# Represents a complete critical box intersection profile
struct BoxIntersections{T}
    bottom::BoxEdge{T}
    top::BoxEdge{T}
    left::BoxEdge{T}
    right::BoxEdge{T}
end

# Tracks the vertex IDs attached to different parts of a critical point box
Base.@kwdef struct CriticalNodeConnections
    left::Vector{Int} = Int[]
    bottom::Vector{Int} = Int[]
    center::Vector{Int} = Int[]
    top::Vector{Int} = Int[]
    right::Vector{Int} = Int[]
end

# Connectivity graph + control nodes
struct CurveGraph{T}
    vertices::Vector{Tuple{T, T}}
    edges::Vector{Tuple{Int, Int}}
    control_nodes::Dict{Int, Vector{Int}}
end

# ---------- PLOT ----------

struct EdgeGroup
    edges::Vector{Tuple{Vector{Float64}, Vector{Float64}}}
    color::String
    width::Float64
end

struct PointGroup
    vertices::Tuple{Vector{Float64}, Vector{Float64}}  # (x, y)
    color::String
    marker::Symbol
end

struct GraphPlotData
    edge_group::EdgeGroup
    point_groups::Vector{PointGroup}
end