# =========================================================================
# GRAPH DECOMPOSITION
# =========================================================================

"""
    connected_components(G::CurveGraph{T}) where T

Decomposes a graph into its connected components.

# Inputs
- `G::CurveGraph{T}`: The unified curve graph.

# Outputs
- `Vector{CurveGraph{T}}`: A list of disjoint subgraphs. Each subgraph retains 
  the mapped coordinates, edges, and control nodes belonging to its component.
"""
function connected_components(G::CurveGraph{T}) where T
    nv = length(G.vertices)

    # Work with vertices to avoid missing isolated vertices
    adj_list = [Int[] for _ in 1:nv]
    for (u, v) in G.edges
        push!(adj_list[u], v)
        push!(adj_list[v], u)
    end

    visited = falses(nv)
    components = Vector{Vector{Int}}()

    # Iterative DFS
    for i in 1:nv
        if !visited[i]
            comp = Int[]
            stack = [i]
            visited[i] = true
            while !isempty(stack)
                curr = pop!(stack)
                push!(comp, curr)
                for neighbor in adj_list[curr]
                    if !visited[neighbor]
                        visited[neighbor] = true
                        push!(stack, neighbor)
                    end
                end
            end
            push!(components, comp)
        end
    end

    # Mapping global vertex ID -> (component ID, local vertex ID)
    vertex_map = Vector{Tuple{Int, Int}}(undef, nv)
    for (cid, comp) in enumerate(components)
        for (local_id, global_id) in enumerate(comp)
            vertex_map[global_id] = (cid, local_id)
        end
    end

    # Group edges
    comp_edges = [Tuple{Int, Int}[] for _ in 1:length(components)]
    for (u, v) in G.edges
        cid_u, local_u = vertex_map[u]
        cid_v, local_v = vertex_map[v]
        push!(comp_edges[cid_u], (local_u, local_v))
    end

    # Reconstruct subgraphs into CurveGraphs
    subgraphs = CurveGraph{T}[]
    for (cid, comp) in enumerate(components)
        new_verts = G.vertices[comp]
        new_edges = comp_edges[cid]
        new_control = Dict{Int, Vector{Int}}()

        # Keep track of control points falling into this specific component
        for (k, v_list) in G.control_nodes
            local_v_list = Int[]
            for global_v in v_list
                c, local_v = vertex_map[global_v]
                if c == cid
                    push!(local_v_list, local_v)
                end
            end
            if !isempty(local_v_list)
                new_control[k] = local_v_list
            end
        end

        push!(subgraphs, CurveGraph{T}(new_verts, new_edges, new_control))
    end

    return subgraphs
end

"""
    group_by_component(G::CurveGraph)

# Outputs
- `Vector{Dict{Int, Vector{Int}}}`: Returns just the control nodes mapped to each connected component.
"""
function group_by_component(G::CurveGraph)
    subgraphs = connected_components(G)
    return [sg.control_nodes for sg in subgraphs if !isempty(sg.vertices)]
end

"""
    number_connected_components(G::CurveGraph)

# Outputs
- `Int`: The count of isolated subgraphs in the parent structure.
"""
number_connected_components(G::CurveGraph) = length(connected_components(G))


# =========================================================================
# GRAPH MERGING
# =========================================================================

"""
    merge_graphs(graphs::Vector{CurveGraph{T}}) where T

Merges a collection of disjoint graphs into a single graph, mathematically fusing
vertices that share common `control_nodes` mappings.

# Inputs
- `graphs::Vector{CurveGraph{T}}`: A list of graphs. `graphs[i].control_nodes[k]`
  should contain the local vertex indices in graph `i` that connect to graph `k`.

# Outputs
- `CurveGraph{T}`: A unified graph struct.
"""
function merge_graphs(graphs::Vector{CurveGraph{T}}) where T
    isempty(graphs) && return CurveGraph{T}(Tuple{T,T}[], Tuple{Int,Int}[], Dict{Int, Vector{Int}}())
    length(graphs) == 1 && return deepcopy(graphs[1])

    Vtot = copy(graphs[1].vertices)
    Etot = copy(graphs[1].edges)

    # Dictionary to track: (graph_id, local_vertex_id) -> merged_global_id
    global_map = Dict{Tuple{Int, Int}, Int}()
    for i in 1:length(graphs[1].vertices)
        global_map[(1, i)] = i
    end

    for i in 2:length(graphs)
        G_i = graphs[i]

        # Step 1: Pre-map common control nodes connected to previously merged graphs
        for (k, common_indices_in_i) in G_i.control_nodes
            if k < i && haskey(graphs[k].control_nodes, i)
                common_indices_in_k = graphs[k].control_nodes[i]
                # Zip ensures we pair corresponding control nodes precisely
                for (idx_i, idx_k) in zip(common_indices_in_i, common_indices_in_k)
                    global_map[(i, idx_i)] = global_map[(k, idx_k)]
                end
            end
        end

        # Step 2: Add completely new (unmapped) vertices from G_i to Vtot
        for idx_i in 1:length(G_i.vertices)
            if !haskey(global_map, (i, idx_i))
                push!(Vtot, G_i.vertices[idx_i])
                global_map[(i, idx_i)] = length(Vtot)
            end
        end

        # Step 3: Append the edges using the newly mapped global indices
        for (u, v) in G_i.edges
            push!(Etot, (global_map[(i, u)], global_map[(i, v)]))
        end
    end

    # Step 4: Preserve external control nodes pointing OUTSIDE the merged structure
    merged_control = Dict{Int, Vector{Int}}()
    num_graphs = length(graphs)

    for i in 1:num_graphs
        for (k, v_list) in graphs[i].control_nodes
            # If `k` points to a graph outside of the ones we just merged
            if k > num_graphs || k <= 0
                mapped_v = [global_map[(i, v)] for v in v_list]
                append!(get!(merged_control, k, Int[]), mapped_v)
            end
        end
    end

    return CurveGraph{T}(Vtot, Etot, merged_control)
end