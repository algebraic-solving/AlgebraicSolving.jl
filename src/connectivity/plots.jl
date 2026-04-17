"""
    distinguishable_colors(n; s=0.8, v=0.9)

Return `n` visually distinct HEX colors using manual HSV->RGB conversion.
No external libraries used.
"""
function distinguishable_colors(n::Int; s::Float64=0.6, v::Float64=0.9)
    n <= 0 && return String[]

    colors = Vector{String}(undef, n)

    for i in 0:n-1
        h = i / n  # [0,1)

        c = v * s
        h6 = h * 6
        x = c * (1 - abs(mod(h6, 2) - 1))
        m = v - c

        r1 = g1 = b1 = 0.0

        if 0 ≤ h6 < 1
            r1, g1, b1 = c, x, 0
        elseif 1 ≤ h6 < 2
            r1, g1, b1 = x, c, 0
        elseif 2 ≤ h6 < 3
            r1, g1, b1 = 0, c, x
        elseif 3 ≤ h6 < 4
            r1, g1, b1 = 0, x, c
        elseif 4 ≤ h6 < 5
            r1, g1, b1 = x, 0, c
        else
            r1, g1, b1 = c, 0, x
        end

        r = Int(round((r1 + m) * 255))
        g = Int(round((g1 + m) * 255))
        b = Int(round((b1 + m) * 255))

        colors[i+1] = @sprintf("#%02X%02X%02X", r, g, b)
    end

    return colors
end


# ---------- CORE BUILDERS ----------

@doc Markdown.doc"""
    build_graph_data(G::CurveGraph{T}; width=3.0, vert=true, color="#FF0000") where T <: Union{Float64,QQFieldElem}

Construct plotting data from a `CurveGraph`.

This function converts the piecewise linear structure encoded in `G` into a
`GraphPlotData` object, separating edges and point groups (vertices and control points),
ready for visualization with plotting libraries such as `Plots.jl` (see documentation
of `GraphPlotData`)

# Arguments
- `G::CurveGraph{T}`: Input graph containing vertices, control nodes, and edges.
- `width::Real=3.0`: Line width used for rendering edges.
- `vert::Bool=true`: If `true`, include vertices in the graph.
- `color::AbstractString="#FF0000"`: Color used for both edges and points.

# Behavior
- Edges are converted into line segments between vertex coordinates.
- Vertices (if enabled) are displayed with marker `:x`.
- Control nodes are grouped and displayed with marker `:+`.

# Returns
- `GraphPlotData`: A structure containing:
  - one `EdgeGroup` for all edges
  - multiple `PointGroup`s for vertices and control nodes

# Notes
- All coordinates are converted to `Float64` for compatibility with plotting backends.
"""
function build_graph_data(G::CurveGraph{T}; width=3.0, vert=true, color="#FF0000")  where T <: Union{Float64,QQFieldElem}
    V, Vcon, E = G.vertices, G.control_nodes, G.edges

    # edges
    edges = [
        begin
            v1, v2 = [map(Float64, V[idx]) for idx in e]
            ([v1[1], v2[1]], [v1[2], v2[2]])
        end
        for e in E
    ]
    edge_group = EdgeGroup(edges, color, width)

    #points
    point_groups = PointGroup[]

    if vert
        vx = map(v -> Float64(v[1]), V)
        vy = map(v -> Float64(v[2]), V)
        push!(point_groups, PointGroup((vx, vy), color, :x))
    end

    for group in Vcon
        hx = map(j -> Float64(V[j][1]), group)
        hy = map(j -> Float64(V[j][2]), group)
        push!(point_groups, PointGroup((hx, hy), color, :+))
    end

    return GraphPlotData(edge_group, point_groups)
end

"""
    build_graph_data(CG::Vector{CurveGraph{T}}; width=3.0, vert=true) where T <: Union{Float64,QQFieldElem}

Construct plotting data for multiple `CurveGraph`s with automatically assigned colors.

# Arguments
- `CG::Vector{CurveGraph{T}}`: Collection of curve graphs to visualize.
- `width::Real=3.0`: Line width used for rendering edges.
- `vert::Bool=true`: If `true`, include vertices in each graph.

# Behavior
- Colors are generated using `distinguishable_colors` to ensure visual contrast.
- Each graph is assigned a unique color consistently across its edges and points.

# Returns
- `Vector{GraphPlotData}`: One plotting data object per input graph.

# Notes
- The number of colors scales with `length(CG)`.
"""
function build_graph_data(CG::Vector{CurveGraph{T}}; width=3.0, vert=true) where T <: Union{Float64,QQFieldElem}
    c = distinguishable_colors(length(CG))
    return [build_graph_data(G; width=width, vert=vert, color=c[i])
            for (i, G) in enumerate(CG)]
end