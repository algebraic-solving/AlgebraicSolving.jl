
export compute_graph, connected_components, number_connected_components, group_by_component, merge_graphs,
 plot_graph, plot_graphs, plot_graph_comp, Bresultant, param_crit_split, build_graphs_data, build_graph_data

 # DEBUG
 export interp_subresultants, mmod_subresultants, subresultants, diff, diff_list, trimat_rand, fact_gcd, isolate_eval, isolate,
 rat_to_arb, evaluate_arb, evaluate_arb_rat, int_coeffs, array_to_poly, parray_asvar, poly_to_array, homogenize, rem_var,
 intersect_biv, num_biv_rat_mod, parray_asvarcoeff, mmod_param_crit, MPolyBuild, change_ringvar

include("datastruct.jl")
include("tools.jl")
include("subresultants-bis.jl")
include("isolateboxes.jl")
include("graph.jl")
include("plots.jl")
include("arbtools.jl")
include("buildpoly.jl")

# =========================================================================
# MULTIPLE DISPATCH WRAPPERS
# =========================================================================

# Base case: No 'C' (control points) provided
compute_graph(f::P, g::P; kwargs...) where {P <: MPolyRingElem} =
    _compute_graph_core(f, g, Vector{Vector{P}}(); kwargs...)

# Case 1: C is a Vector of Vectors
compute_graph(f::P, g::P, C::Vector{Vector{P}}; kwargs...) where {P <: MPolyRingElem} =
    _compute_graph_core(f, g, C; kwargs...)

# Case 2: C is a Dictionary
function compute_graph(f::P, g::P, C::Dict{Int, Vector{P}}; kwargs...) where {P <: MPolyRingElem}
    graph = _compute_graph_core(f, g, collect(values(C)); kwargs...)
    # Remap Vcon dictionary keys to match original dict keys
    mapped_vcon = Dict{Int, Vector{Int}}(k => graph.control_nodes[i] for (i, k) in enumerate(keys(C)))
    return CurveGraph(graph.vertices, graph.edges, mapped_vcon)
end

# Case 3: C is a single Vector
compute_graph(f::P, g::P, C::Vector{P}; kwargs...) where {P <: MPolyRingElem} =
    _compute_graph_core(f, g, [C]; kwargs...)

# =========================================================================
# CORE IMPLEMENTATION
# =========================================================================

# Input:
# The space curve is given by f(x,y) = 0 and (df/dy)(x,y)*z = g(x,y)
# it is assumed in generic position of put it using generic parameter
# for now f is assumed to be square-free (future feature)
# C is a list of parametrization [p,a,b] of plane pts i.e. such that p(x)=0,y=a(x)/b(x)
# these points must be on the curve
# Output: Computes a graph homeomorphic to the curve, identitying point in C
function _compute_graph_core(f::P, g::P, C::Vector{Vector{P}};
                             generic=true, precx=150, v=0, int_coeff=true, outf=true) where {P <: MPolyRingElem}

    R = parent(f)
    x, y = gens(R)
    intC = int_coeff ? int_coeffs : identity

    # Pre-processing the input
    f, g = intC([f, g])
    changemat = generic ? [1 0; 0 1] : trimat_rand(QQ, 2, range=-100:100)
    f = evaluate(f, collect(changemat * [x; y]))
    precx = max(2, precx)
    v > 1 && println(f)

    v > 0 && println("Compute parametrization of critical pts...")
    @iftime (v > 0) params = param_crit_split(f, g, v=v-1, detect_app=true)
    for i in 1:length(C)
        params[-i] = [ [intC(C[i][1])], C[i][2], C[i][3] ]
    end

    v > 0 && println("Computing insulating critical boxes")
    @iftime (v > 0) LBcrit, Lprecx = insulate_crit_boxes(f, params, precx, v=v-1)

    v > 0 && println("Compute intersections with critical boxes..")
    @iftime (v > 0) LPCside, LnPCside = intersect_vertical_boxes(f, params, LBcrit, Lprecx, outf, v=v-1)

    # Critical values and their order
    xcrit = Dict(i => [LBcrit[i][j][1] for j in eachindex(LBcrit[i])] for i in keys(LBcrit))
    xcritpermut = order_permut2d(xcrit)

    # Graph data
    Vert = Tuple{QQFieldElem, QQFieldElem}[] # List of points (x,y)
    Edg = Tuple{Int, Int}[] # List of tuples (idx, idy)
    Vcon = Dict{Int, Vector{Int}}(k => Int[] for k in keys(C)) # Index of control vertices

    Lapp_isolated = Tuple{Int, Int}[]
    Lapp_nodes = Tuple{Int, Int}[]

    # Keep track of connection between all computed points when increasing x
    Corr = Dict(m => [CriticalNodeConnections() for _ in xcrit[m]] for m in keys(xcrit))

    for (ind, (i, j)) in enumerate(xcritpermut)
        i1, j1 = ind > 1 ? xcritpermut[ind - 1] : (0, 0) # previous
        i2, j2 = ind < length(xcritpermut) ? xcritpermut[ind + 1] : (0, 0) # next

        # Current cylinder above xcrit[i][j]
        PCside = LPCside[i][j]
        I2L_points = ind < length(xcritpermut) ? LPCside[i2][j2].left.points : nothing

        # Midpoints of the critical box itself
        xcmid = sum(LBcrit[i][j][1]) / 2
        ycmid = sum(LBcrit[i][j][2]) / 2

        nI_left = PCside.left.indices_inside
        nI_right = PCside.right.indices_inside
        ymincrit = minimum(vcat(nI_left, nI_right, [length(PCside.left.points) + 1]))

        # 1. Left side vertices
        # one-to-one connection with previous right side
        if ind > 1
            for k in 1:length(PCside.left.points)
                push!(Corr[i][j].left, Corr[i1][j1].right[k])
            end
        else
            for k in 1:length(PCside.left.points)
                xval = xcrit[i][j][1]
                yval = sum(PCside.left.points[k]) / 2
                push!(Vert, (xval, yval))
                push!(Corr[i][j].left, length(Vert))
            end
        end

        # 2. Right side vertices
        if ind < length(xcritpermut)
            for k in 1:length(PCside.right.points)
                xval = (xcrit[i][j][2] + xcrit[i2][j2][1]) / 2
                yval = sum(PCside.right.points[k] + I2L_points[k]) / 4
                push!(Vert, (xval, yval))
                push!(Corr[i][j].right, length(Vert))
            end
        else
            for k in 1:length(PCside.right.points)
                xval = xcrit[i][j][2]
                yval = sum(PCside.right.points[k]) / 2
                push!(Vert, (xval, yval))
                push!(Corr[i][j].right, length(Vert))
            end
        end

        # 3. Below the critical point
        # TODO: many vertices and edges could be removed here
        for k in 1:ymincrit-1
            yval = sum(PCside.left.points[k] + PCside.right.points[k]) / 4
            push!(Vert, (xcmid, yval))
            push!(Corr[i][j].bottom, length(Vert))
            push!(Edg, (Corr[i][j].left[k], length(Vert)))
            push!(Edg, (length(Vert), Corr[i][j].right[k]))
        end

        # 4. The critical point
        if i == 0 && isempty(nI_left) && isempty(nI_right)
            push!(Lapp_isolated, (i, j))
        elseif i == 0
            push!(Edg, (Corr[i][j].left[nI_left[1]], Corr[i][j].right[nI_right[2]]))
            push!(Edg, (Corr[i][j].left[nI_left[2]], Corr[i][j].right[nI_right[1]]))
            push!(Lapp_nodes, (i, j))
        else
            push!(Vert, (xcmid, ycmid))
            push!(Corr[i][j].center, length(Vert))

            for k in nI_left; push!(Edg, (Corr[i][j].left[k], length(Vert))); end
            for k in nI_right; push!(Edg, (length(Vert), Corr[i][j].right[k])); end

            if i < 0 # Control point
                push!(Vcon[-i], length(Vert))
            end
        end

        # 5. Above the critical point
        for k = (length(PCside.left.points) - length(nI_left) - ymincrit + 1):-1:1
            yval = sum(PCside.left.points[end-k+1] + PCside.right.points[end-k+1]) / 4
            push!(Vert, (xcmid, yval))
            push!(Corr[i][j].top, length(Vert))
            push!(Edg, (Corr[i][j].left[end-k+1], length(Vert)))
            push!(Edg, (length(Vert), Corr[i][j].right[end-k+1]))
        end
    end

    # Post processing the data
    if !generic
        for i in eachindex(Vert)
            vx, vy = Vert[i]
            Vert[i] = Tuple(changemat * [vx, vy])
        end
    end

    if outf
        Vert = [ (Float64(x), Float64(y)) for (x,y) in Vert ]
        return CurveGraph{Float64}(Vert, Edg, Vcon)
    end

    # Return the unified data structure
    return CurveGraph{QQFieldElem}(Vert, Edg, Vcon)
end