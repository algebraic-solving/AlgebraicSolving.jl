
export compute_graph, connected_components, number_connected_components, group_by_component, merge_graphs,
 plot_graph, plot_graphs, plot_graph_comp, Bresultant, param_crit_split

 # DEBUG
 export interp_subresultants, mmod_subresultants, subresultants, diff, diff_list, trimat_rand, fact_gcd, isolate_eval, isolate,
 rat_to_arb, evaluate_arb, evaluate_arb_rat, int_coeffs, array_to_poly, parray_asvar, poly_to_array, homogenize, rem_var,
 intersect_biv, num_biv_rat_mod, parray_asvarcoeff, mmod_param_crit

include("tools.jl")
include("subresultants.jl")
include("isolateboxes.jl")
include("graph.jl")
include("plots.jl")
include("arbtools.jl")
include("buildpoly.jl")

function compute_graph(f::P, g::P, C::Vector{Vector{P}}=Vector{Vector{P}}(); generic=true, precx = 150, v=0, int_coeff=true, outf=true)  where (P <: MPolyRingElem)
    R = parent(f)
    x, y = gens(R)
    intC = int_coeff ? int_coeffs : identity
    # Pre-processing the input
    f, g  = intC([f,g])
    changemat = generic ? [1 0; 0 1] : trimat_rand(QQ, 2, range=-100:100)
    f = evaluate(f, collect(changemat*[x; y]));
    v > 1 && println(f)
    precx = precx <= 1 ? 2 : precx

    v > 0 && println("Compute parametrization of critical pts...")
    @iftime (v > 0) params = param_crit_split(f,g,v=v-1, detect_app=true)
    #@iftime (v > 0) params = mmod_param_crit(f, g,v=v-1, detect_app=true)

    for i in 1:length(C)
        params[-i] = [ [C[i][1] |> intC], C[i][2], C[i][3] ]
    end

    v > 0 && println("Computing isolating critical boxes")
    LBcrit = Dict()
    @iftime (v > 0) xcritpermut, precx = compute_crit_boxes(params, LBcrit, precx, v)
    xcrit = Dict(i => [LBcrit[i][j][1] for j in eachindex(LBcrit[i])] for i in eachindex(LBcrit))

    # println("LBcrit: ", LBcrit)

    v > 0 && println("\nTest for identifying singular boxes")
@iftime (v > 0) begin
    ########################################################
    # For each mult of sing pts (keys) give the corresponding factors of the resultant (values)
    #########################################################
    arbField = ArbField(precx)
    Lfyk = diff_list(f, 2, max(maximum(eachindex(LBcrit)),2))
    for ind in eachindex(LBcrit)
        ind < 0 &&  continue # control pts
        m = ind <= 1 ? 2 : ind # nodes and extrems have mult 2
        noisolate = false
        for _ in 1:5
            noisolate = false
            for j in eachindex(LBcrit[ind])
                pcrit = [ rat_to_arb(c, arbField) for c in LBcrit[ind][j] ]
                # Check if the the mult(pcrit)-th derivative of f vanishes on pcrit
                if contains_zero(evaluate(Lfyk[m+1], pcrit))
                    # TODO refine only boxes associated to multiplicity m
                    precx *= 2
                    (v > 0) && println("Refine singular boxes of index $ind to $precx bits precision")
                    xcritpermut, precx = compute_crit_boxes(params, LBcrit, precx, v-1)
                    xcrit = Dict(i => [LBcrit[i][j][1] for j in eachindex(LBcrit[i])] for i in eachindex(LBcrit))
                    noisolate = true
                    break
                end
            end
            noisolate || break
        end
        if noisolate
            error("Problem in singular boxes refinement")
        end
    end
end

    v > 0 && println("\nCompute intersections with critical boxes..")
@iftime (v > 0) begin
    ## TODO : Refine only the intervals that need to be refined? only the factor that defines the intervals
    LPCside = Dict{Int,Any}()
    ndig = maximum([ndigits(length(LBcrit[i])) for i in eachindex(LBcrit)])
    for i in eachindex(LBcrit)
        LPCside[i] = Array{Any}(undef, length(LBcrit[i]))
        # Data printing part
        ndigi = ndigits(length(LBcrit[i]))
        Ptype = (i > length(LBcrit)-length(C)) ? "Pcon" : "mult"
        ###
        precxtmp = precx
        compt = 0
        while compt < 5
            flag = false
            # println(LBcrit[i])
            for j in eachindex(LBcrit[i])
                v > 0 && print("$Ptype=$i ; $(j)/$(length(LBcrit[i]))$(repeat(" ", ndig-ndigi+1))pts","\r")
                pcside = intersect_box(f, LBcrit[i][j], prec=precxtmp,v=v)
                npcside = [length(n) for (I, n) in pcside]
                # println("($i, $j): $npcside")
                if (i == 1 && sum(npcside) > 2) || (i != 1 && sum(npcside[1:2]) != 0)
                    precxtmp *= 2
                    v > 0 && println("\nRefine boxes along x-axis to $precxtmp bits precision")
                    refine_xboxes(params[i][1], LBcrit[i], precxtmp)
                    flag = true
                    compt += 1
                    break
                end
                LPCside[i][j] = pcside
            end
            flag || break
        end
        v > 0 && println("")
        if compt >= 5
            error("Problem in computing intersections with boxes")
        end
    end
    LnPCside = Dict(i => [[length(indI) for (L, indI) in PB] for PB in LPCside[i]] for i in eachindex(LPCside))

    # Update extreme boxes
    if haskey(LBcrit, 1)
        for j in eachindex(LBcrit[1])
            # If the curve does not intersect the box only on vertical sides
            if !(LnPCside[1][j][1:2] == [0, 0])
                PCside, nPCside = LPCside[1][j], LnPCside[1][j]
                I = [ l[1] for l in PCside[3:end] ]
                nI = [ l[2] for l in PCside[3:end] ]
                # Locate the orientation of the extreme point
                # s is the index on the side where there are more branches
                # s=1: left; s=2: right
                s = argmax([length(I[1]), length(I[2])])
                # Ordinate range of the extreme point
                ycrit = LBcrit[1][j][2]
                # If it intersects on the bottom side
                if nPCside[1] == 1
                    # yinf: the intersection with the vertical side just below the extreme point
                    yinf = maximum([i for (i, yy) in pairs(I[s]) if yy[1] < ycrit[1]])
                    # We vertically enlarge the box until it intersects on the horizontal side
                    push!(LPCside[1][j][s + 2][2], yinf)
                    LPCside[1][j][1][2] = []
                    # We update the intersection numbers
                    LnPCside[1][j][s + 2] += 1
                    LnPCside[1][j][1] = 0
                end
                # If it intersects on the top side
                if nPCside[2] == 1
                    # ymax: the intersection with the vertical side just above the extreme point
                    ymax = minimum([i for (i, yy) in pairs(I[s]) if yy[2] > ycrit[2]])
                    # We vertically enlarge the box until it intersects on the horizontal side
                    push!(LPCside[1][j][s + 2][2], ymax)
                    LPCside[1][j][2][2] = []
                    # We update the intersection numbers
                    LnPCside[1][j][s + 2] += 1
                    LnPCside[1][j][2] = 0
                end
            end
        end
    end
end

    # Graph computation
    fct, typeout = outf ? (Float64, Float64) : (identity, QQFieldElem)

    # List of vertex: at index i is the tuple of the coordinates of the i-th vertex
    Vert = Vector{Tuple{typeout, typeout}}()
    # List of edges: a tuple (i,j) is an edge between the i-th and j-th vertices
    Edg = Vector{Tuple{Int, Int}}()
    Corr = Dict( m => [[[], [[], [], []], []] for j in xcrit[m] ] for m in eachindex(xcrit))
    Vcon = [ [] for _ in 1:length(C) ]

    # Keep track of processed isolated & apparent singularities
    Lapp = [[],[]]

    for ind in eachindex(xcritpermut)
        i, j = xcritpermut[ind]

        if ind > 1
            i1, j1 = xcritpermut[ind - 1]
        end
        if ind < length(xcritpermut)
            i2, j2 = xcritpermut[ind + 1]
            I2L, = LPCside[i2][j2][3]
        end

        PCside = LPCside[i][j]
        I = [ l[1] for l in PCside[3:end] ]
        nI = [ l[2] for l in PCside[3:end] ]

        xcmid = sum(LBcrit[i][j][1])//2 |> fct
        ycmid = sum(LBcrit[i][j][2])//2 |> fct

        ymincrit = minimum(vcat(nI[1], nI[2], [length(I[1])+1]))
        # Construct vertices
        ###########################
        # On the vertical left side
        if ind > 1
            for k in 1:length(I[1])
                push!(Corr[i][j][1], Corr[i1][j1][3][k])
            end
        else
            for k in 1:length(I[1])
                push!(Vert, (xcrit[i][j][1], sum(I[1][k])//2) .|> fct)
                push!(Corr[i][j][1], length(Vert))
            end
        end
        ###########################
        # On the vertical right side
        if ind < length(xcritpermut)
            for k in 1:length(I[2])
                push!(Vert, ((xcrit[i][j][2] + xcrit[i2][j2][1])//2, sum(I[2][k] + I2L[k])//4) .|> fct)
                push!(Corr[i][j][3], length(Vert))
            end
        else
            for k in 1:length(I[2])
                push!(Vert, (xcrit[i][j][2], sum(I[2][k])//2) .|> fct)
                push!(Corr[i][j][3], length(Vert))
            end
        end
        ###########################
        # Below the critical point
        for k in 1:ymincrit-1
            push!(Vert, (xcmid, sum(I[1][k] + I[2][k])// 4) .|> fct)
            push!(Corr[i][j][2][1], length(Vert))
            push!(Edg, (Corr[i][j][1][k], length(Vert)))  # left
            push!(Edg, (length(Vert), Corr[i][j][3][k]))  # right
        end
        ###########################
        # The critical point
        ##########################
        # If we are dealing with a isolated node
        if i in [0,2] && isempty(nI[1]) && isempty(nI[2])
                push!(Lapp[1], (i,j))
        # If we are dealing with an apparent singularity
        elseif i == 0
            # We connect the pairwise opposite branches nI[1][1][i] and nI[1][2][i+1 mod 2], i=1,2
            push!(Edg, (Corr[i][j][1][nI[1][1]], Corr[i][j][3][nI[2][2]]))
            push!(Edg, (Corr[i][j][1][nI[1][2]], Corr[i][j][3][nI[2][1]]))
            push!(Lapp[2], (i,j))
        else
            # We can add the vertex
            push!(Vert, (xcmid, ycmid))
            push!(Corr[i][j][2][2], length(Vert))
            # We connect to the vertical left side of the critical box
            for k in nI[1]
                push!(Edg, (Corr[i][j][1][k], length(Vert)))
            end
            # We connect to the vertical right side of the critical box
            for k in nI[2]
                push!(Edg, (length(Vert), Corr[i][j][3][k]))
            end
            if i < 0
                # If this is a control point
                push!(Vcon[-i], length(Vert))
            end
        end
        ###########################
        # Above the critical point
        for k=(length(I[1]) - length(nI[1]) - ymincrit+1):-1:1
            push!(Vert, (xcmid, sum(I[1][end - k + 1] + I[2][end - k + 1])//4) .|> fct)
            push!(Corr[i][j][2][3], length(Vert))
            push!(Edg, (Corr[i][j][1][end - k + 1], length(Vert)))  # left
            push!(Edg, (length(Vert), Corr[i][j][3][end - k + 1]))  # right
        end
    end

    if v > 0 && Lapp != [[],[]]
        println("Removed apparent singularities: $(length(Lapp[1])) isolated  $(length(Lapp[2])) nodes")
    end

    # Operate inverse change of variable if necessary
    if !(generic)
        for i in eachindex(Vert)
            vx,vy = Vert[i]
            vx,vy = fct.(changemat*[vx, vy])
            Vert[i] = (vx, vy)
        end
    end

    if length(C)>0
        return (Vert, Edg), Vcon
    else
        return Vert, Edg
    end
end

function compute_graph(f::P, g::P,  C::Dict{Int,Vector{P}}; generic=true, precx = 150, v=0, arb=true, int_coeff=false, outf = true) where (P <: MPolyRingElem)
    G, Vcon = compute_graph(f, g, collect(values(C)), generic=generic, precx=precx, v=v, arb=arb, int_coeff=int_coeff, outf=outf)
    VC = Dict{Int64, Vector{Int}}([ k => Vcon[i] for (i,k) in enumerate(keys(C)) ])
    return G, VC
end

function compute_graph(f::P, g::P, C::Vector{P}; generic=true, precx = 150, v=0, arb=true, int_coeff=false, outf = true) where (P <: MPolyRingElem)
    G, (Vcon,) = compute_graph(f, g, [C], generic=generic, precx=precx, v=v, arb=arb, int_coeff=int_coeff, outf = outf)
    return G, Vcon
end

################################################################
function group_by_component(f::P, C::Vector{P}; generic=true, precx=150, v=0, arb=true, int_coeff=false) where {P <: MPolyRingElem}
    # Compute a graph homeomorphic to Z(F) and return the vertices identified by Z(C)
    G, Vcon = compute_graph(f, C, generic=generic, precx=precx, v=v, arb=arb, int_coeff=int_coeff)
    # Compute the partition of the vertices in Vcon according to the connected component in G
    sort!(Vcon)
    CVcon = group_by_component(G, Vcon)

    # Convert the partition into abscissa order
    index_map = Dict((val, idx) for (idx, val) in enumerate(Vcon))
    return [map(v -> index_map[v], C) for C in CVcon]
end

function plot_graph(f::P, g::P, C::Vector{Vector{P}}=Vector{Vector{P}}(); generic=true, precx = 150, v=0, arb=true, int_coeff=true, outf=true)  where (P <: MPolyRingElem)
    plot_graph(compute_graph(f, g, C, generic=generic, precx=precx, v=v, arb=arb, int_coeff=int_coeff, outf=outf))
end

function plot_graph(f::P, g::P, C; generic=true, precx = 150, v=0, arb=true, int_coeff=true, outf=true)  where (P <: MPolyRingElem)
    plot_graph(compute_graph(f, g, C, generic=generic, precx=precx, v=v, arb=arb, int_coeff=int_coeff, outf=outf))
end

function plot_graph_comp(f::P, g::P, C::Vector{Vector{P}}=Vector{Vector{P}}(); generic=true, precx = 150, v=0, arb=true, int_coeff=true, outf=true)  where (P <: MPolyRingElem)
    plot_graph_comp(compute_graph(f, g, C, generic=generic, precx=precx, v=v, arb=arb, int_coeff=int_coeff, outf=outf))
end

function plot_graph_comp(f::P, g::P, C; generic=true, precx = 150, v=0, arb=true, int_coeff=true, outf=true)  where (P <: MPolyRingElem)
    plot_graph_comp(compute_graph(f, g, C, generic=generic, precx=precx, v=v, arb=arb, int_coeff=int_coeff, outf=outf))
end

function number_connected_components(f::P, g::P, C::Vector{Vector{P}}=Vector{Vector{P}}(); generic=true, precx = 150, v=0, arb=true, int_coeff=true, outf=true)  where (P <: MPolyRingElem)
    number_connected_components(compute_graph(f, g, C, generic=generic, precx=precx, v=v, arb=arb, int_coeff=int_coeff, outf=outf))
end

function number_connected_components(f::P, g::P, C; generic=true, precx = 150, v=0, arb=true, int_coeff=true, outf=true)  where (P <: MPolyRingElem)
    number_connected_components(compute_graph(f, g, C, generic=generic, precx=precx, v=v, arb=arb, int_coeff=int_coeff, outf=outf))
end