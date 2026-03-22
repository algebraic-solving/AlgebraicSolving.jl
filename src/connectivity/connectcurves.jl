#println("\nLoading libraries and data..")
#using Nemo
#using Plots, Colors
#pythonplot()

export compute_graph, connected_components, number_connected_components, group_by_component, merge_graphs,
 plot_graph, plot_graphs, plot_graph_comp, Bresultant, param_crit_split

 # DEBUG
 export interp_subresultants, mmod_subresultants, subresultants, diff, diff_list, trimat_rand, fact_gcd, isolate_eval, isolate,
 rat_to_Arb, evaluate_Arb, evaluate_Arb_rat, int_coeffs, array_to_poly, parray_asvar, poly_to_array, homogenize, rem_var,
 intersect_biv, num_biv_rat_mod, parray_asvarcoeff, mmod_param_crit

include("tools.jl")
include("subresultants.jl")
include("isolate.jl")
include("boxes.jl")
include("graph.jl")
include("plots.jl")
include("arbtools.jl")
include("src/resultant/bresultant.jl")
include("buildpoly.jl")

function compute_graph(f::P, g::P, C::Vector{Vector{P}}=Vector{Vector{P}}(); generic=true, precx = 150, v=0, arb=true, int_coeff=true, outf=true)  where (P <: MPolyRingElem)
    R = parent(f)
    x, y = gens(R)
    intC = int_coeff ? int_coeffs : identity
    # Pre-processing the input
    f, g  = intC([f,g])
    changemat = generic ? [1 0; 0 1] : trimat_rand(QQ, 2, range=-100:100)
    f = evaluate(f, collect(changemat*[x; y]));
    v > 1 && println(f)

    v > 0 && println("Compute parametrization of critical pts...")
    @iftime (v > 0) params = param_crit_split(f,g,v=v-1, detect_app=true)
    #@iftime (v > 0) params = mmod_param_crit(f, g,v=v-1, detect_app=false)
    #=params = Bresultant(f, derivative(f,2), bspath="AlgebraicSolving.jl/src/connectivity/src/resultant/bresultant",v=2);
    params = change_ringvar.(params, Ref([:x,:y]))
    params = [ [[p[1]], p[3], p[2]] for p in params]
    params = Dict(1=>params[1], -1=>params[2])=#

    for i in 1:length(C)
        params[-i] = [ [C[i][1] |> intC], C[i][2], C[i][3] ]
    end
    #println(keys(params))
    if arb
        compt = 0
        while compt < 5
            try
                # TODO : check that no overlap between different isolations
                v > 0 && println("\nIsolating critical values with precision ", precx,"..")
            @iftime (v > 0) begin
                xcrit = Dict(p[1]=> reduce(vcat, [isolate(pp, prec=precx, software="msolve") for pp in p[2][1]]) for p in params)
                # Enlarge exact isolating root intervals
                for i in eachindex(xcrit)
                    for j in eachindex(xcrit[i])
                        if xcrit[i][j][1]==xcrit[i][j][2]
                            xcrit[i][j] = [xcrit[i][j][1]-1//ZZ(1)<<precx, xcrit[i][j][1]+1//ZZ(1)<<precx]
                        end
                    end
                end
                xcritpermut = order_permut2d(xcrit);
            end

                v > 0 && println("\nComputing isolating critical boxes using Arb with precision ",max(precx,150),"..")
            @iftime (v > 0) begin
                precArb = precx
                Pcrit = Dict( i => [[xc, evaluate_Arb(params[i][2], params[i][3], rat_to_Arb(xc, precArb))] for xc in xcrit[i]] for i in eachindex(xcrit))
                LBcrit =Dict( i=> [[ map(QQ, pc[1]), map(QQ, Arb_to_rat(pc[2])) ]  for pc in pcrit] for (i, pcrit) in Pcrit)
                # Enlarge exact isolating box vertical side (horizontal lines)
                for i in eachindex(LBcrit)
                    for j in eachindex(LBcrit[i])
                        if LBcrit[i][j][2][1]==LBcrit[i][j][2][2]
                            LBcrit[i][j][2] = [LBcrit[i][j][2][1]-1//ZZ(2)<<precx, LBcrit[i][j][2][1]+1//ZZ(2)<<precx]
                        end
                    end
                end
            end
            catch
                precx *= 2
                v > 0 && println("\nRefine x-precision to $precx")
                compt += 1
            else
                break
            end
        end
        if compt >= 5
            error("Problem in isolating critical boxes")
        end
    else
        v > 0 && println("\nCompute critical boxes with msolve with precision ", precx,"..")
    @iftime (v > 0) begin
        # TODO: parameter -I
        LBcrit = Dict(p[1]=>reduce(vcat,[ sort(real_solutions(AlgebraicSolving.Ideal([pp,  p[2][3]*y-p[2][2]]), precision=precx,interval=true),by=t->t[1])  for pp in p[2][1] ]) for p in params)
        # Enlarge each exact isolating box side
        for i in eachindex(LBcrit)
            for j in eachindex(LBcrit[i])
                for k in eachindex(LBcrit[i][j])
                    if LBcrit[i][j][k][1]==LBcrit[i][j][k][2]
                        LBcrit[i][j][k] = [LBcrit[i][j][k][1]-1//ZZ(2-k)<<precx, LBcrit[i][j][k][1]+1//ZZ(2-k)<<precx]
                    end
                end
            end
        end
        xcrit = Dict(lbcrit[1]=>[ B[1] for B in lbcrit[2] ] for lbcrit in LBcrit)
        xcritpermut = order_permut2d(xcrit);
    end
    end
    #println(xcrit)

    v > 0 && println("\nTest for identifying singular boxes")
@iftime (v > 0) begin
    ########################################################
    # For each mult of sing pts (keys) give the corresponding factors of the resultant (values)
    #########################################################
    Lfyk = diff_list(f, 2, max(maximum(eachindex(LBcrit)),2))
    for ind in eachindex(LBcrit)
        ind>1 || ind==0 ||  continue # extreme and control pts
        m = ind==0 ? 2 : ind # nodes have mult 2
        while true
            flag = false
            for j in eachindex(LBcrit[ind])
                pcrit = [ rat_to_Arb(c, precx) for c in LBcrit[ind][j] ]
                # Check if the the mult(pcrit)-th derivative of f vanishes on pcrit
                if contains_zero(evaluate(Lfyk[m+1], pcrit))
                    (v > 0) && println("Refine singular boxes of multiplicity ", m)
                    # TODO refine_boxes(params[k], LBcrit[k])
                    #flag = true
                    break
                end
            end
            flag || break
        end
    end
end

    v > 0 && println("\nCompute intersections with critical boxes..")
@iftime (v > 0) begin
    # Could be improved by handling nodes (or even any ordinary sing) as extreme boxes:
    # when npcside = [2,2,0,0] just take nearest below and above
    # intersections b with the curves on the vertical sides
    # and change into npcside = [0,0,2,2]
    ## TODO : Refine only the intervals that need to be refined
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
            for j in eachindex(LBcrit[i])
                v > 0 && print("$Ptype=$i ; $(j)/$(length(LBcrit[i]))$(repeat(" ", ndig-ndigi+1))pts","\r")
                pcside = intersect_box(f, LBcrit[i][j], prec=precxtmp,v=v)
                npcside = [length(n) for (I, n) in pcside]
                #println("($i, $j): $npcside")
                if (i == 1 && sum(npcside) > 2) || (i != 1 && sum(npcside[1:2]) != 0)
                    precxtmp *= 2
                    v > 0 && println("\nRefine boxes along x-axis to precision ", precxtmp)
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

    #return LnPCside

    #println("\nGraph computation")
    # Would be nice to have only one intermediate fiber (take the average of abscissa and ordinates) for plot
    # And even remove this fiber for the graph
    fct, typeout = outf ? (Float64, Float64) : (identity, QQFieldElem)

    Vert = Vector{Tuple{typeout, typeout}}()
    Edg = Vector{Tuple{Int, Int}}()
    Corr = Dict( m => [[[], [[], [], []], []] for j in xcrit[m] ] for m in eachindex(xcrit))
    Viso = []
    Vcon = [ [] for _ in 1:length(C) ]

    Lapp = [[],[]]

    for ind in eachindex(xcritpermut)
        i, j = xcritpermut[ind]

        if ind > 1
            i1, j1 = xcritpermut[ind - 1]
        end
        if ind < length(xcritpermut)
            i2, j2 = xcritpermut[ind + 1]
            I2L, nI2L = LPCside[i2][j2][3]
        end

        PCside, nPCside = LPCside[i][j], LnPCside[i][j]
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
        #println()
        #println(map(length,I))
        #println(nI[1],", ", nI[2], ", ", [length(I[1])+1])
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
                #pass
                # We can add the isolated  vertex
                # push!(Vert, [xcmid, ycmid])
                # push!(Corr[i][j][2][1], length(Vert))
                # We will subsequently add the vertex in the graph
                # push!(Viso, length(Vert))
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