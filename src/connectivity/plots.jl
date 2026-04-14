
## Plot functions
function plot_graph(G; color="red", width=3, vert=true, subplt=false)
    if !subplt
        plot(legend=false)
    end
    V, Vcon, E = G.vertices, G.control_nodes, G.edges
    #col = distinguishable_colors(length(Vemph)+2)
    for e in E
        v1, v2 = [ map(Float64, V[ee]) for ee in e ]
        plot!([v1[1], v2[1]], [v1[2], v2[2]], lc=color, lw=width)
    end
    if vert
        scatter!( map(Float64, [v[1] for v in V]),  map(Float64, [v[2] for v in V]), mc="black", m=:diamond)
    end
    for i in 1:length(Vcon)
        scatter!( map(Float64, [V[j][1] for j in Vcon[i]]),  map(Float64, [V[j][2] for j in Vcon[i]]), mc=color, m=:diamond)
    end
    subplt || gui()
end

function plot_graphs(CG; width=3, vert=true, subplt=false)
    if !subplt
        plot(legend=false)
    end
    col = distinguishable_colors(length(CG)+2)
    for (j, G) in enumerate(CG)
        plot_graph(G, color=col[j+2], vert=vert, subplt=true)
    end
    subplt || gui()
    #savefig("/home/remi/Documents/gittravail/test.html")
    #replace_width_height_in_file("/home/remi/Documents/gittravail/test.html", 1000, 750)
end

function plot_graph_comp(G; width=3, vert=true, subplt=false)
    plot(legend=false)
    CG = connected_components(G)
    plot_graphs(CG, width=width, vert=vert, subplt=true)
    subplt || gui()#savefig("/tmp/test.html")
end

function replace_width_height_in_file(filename, x, y)
    # Open the file and read its content as a string
    content = read(filename, String)

    # Build the patterns to search for and the replacements
    width_pattern = "\"width\": 600"
    height_pattern = "\"height\": 400"

    # Create the replacement strings with new values of x and y
    new_width = "\"width\": $x"
    new_height = "\"height\": $y"

    # Replace the last occurrence of "width": 600
    last_width_pos = findlast(occursin(width_pattern), Base.split(content, "\n"))
    if last_width_pos !== nothing
        content = replace(content, width_pattern => new_width; count=1)
    end

    # Replace the last occurrence of "height": 400
    last_height_pos = findlast(occursin(height_pattern), Base.split(content, "\n"))
    if last_height_pos !== nothing
        content = replace(content, height_pattern => new_height; count=1)
    end

    # Write the modified content back to the file
    write(filename, content)
end

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