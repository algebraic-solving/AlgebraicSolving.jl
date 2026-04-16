"""
    distinguishable_colors(n; s=0.8, v=0.9)

Return `n` visually distinct HEX colors using manual HSV→RGB conversion.
No external libraries used.
"""
function distinguishable_colors(n::Int; s::Float64=0.8, v::Float64=0.9)
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
"""
# How to plot:
using Plots
function plot_graph(P)
  plot(legend=false)
  E = P.edge_group
  plot!(E.edges, color=E.color, width = E.width)
  [scatter!(V.vertices, color=V.color, marker = V.marker) for V in P.point_groups]
  gui()
end

plot_graph(build_graphs_data(G))
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

# ---------- MULTI GRAPHS ----------

"""
# How to plot:
using Plots
function plot_graphs(CP)
  plot(legend=false)
  for P in CP
      E = P.edge_group
      plot!(E.edges, color=E.color, width = E.width)
      [ scatter!(V.vertices, color=V.color, marker = V.marker) for V in P.point_groups ]
  end
  gui()
end

plot_graphs(build_graph_data(CG))
"""

function build_graph_data(CG::Vector{CurveGraph{T}}; width=3.0, vert=true) where T <: Union{Float64,QQFieldElem}
    c = distinguishable_colors(length(CG))
    return [build_graph_data(G; width=width, vert=vert, color=c[i])
            for (i, G) in enumerate(CG)]
end