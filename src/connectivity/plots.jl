
# ---------- CORE BUILDERS ----------
# How to plot:
# function plot_graph(P)
#   plot(legend=false)
#   E = P.edge_group
#   plot!(E.edges, color=E.color, width = E.width)
#   scatter!(V.vertices, color=V.color, marker = V.marker) for V in P.point_groups ]
#   gui()
# end
#
# plot_graph(build_graphs_data(G))

function build_graph_data(G; width=3.0, vert=true, color=rand(1:666))
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
        push!(point_groups, PointGroup((vx, vy), 17, :diamond))
    end

    for group in Vcon
        hx = map(j -> Float64(V[j][1]), group)
        hy = map(j -> Float64(V[j][2]), group)
        push!(point_groups, PointGroup((hx, hy), color, :diamond))
    end

    return GraphPlotData(edge_group, point_groups)
end

# ---------- MULTI GRAPHS ----------

# How to plot:
# function plot_graphs(P)
#   plot(legend=false)
#   for P in CP
#       E = P.edge_group
#       plot!(E.edges, color=E.color, width = E.width)
#       [ scatter!(V.vertices, color=V.color, marker = V.marker) for V in P.point_groups ]
#   end
#   gui()
# end
#
# plot_graphs(build_graphs_data(CG))

function build_graphs_data(CG; width=3.0, vert=true)
    c = randperm(666)[1:length(CG)] # Distinct numbers
    return [build_graph_data(G; width=width, vert=vert, color=c[i])
            for (i, G) in enumerate(CG)]
end