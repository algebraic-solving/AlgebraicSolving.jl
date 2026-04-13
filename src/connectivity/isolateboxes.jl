const Side = Tuple{Vector{QQFieldElem}, Vector{Int}}
const Sides = NTuple{4, Side}
const Boxes = Vector{Side}


in_inter(I, J) =  J[1] <= I[1] && I[2] <= J[2]
overlap_inter(I,J) = max(I[1], J[1]) <= min(I[2], J[2])

 # Function that isolate roots of a multivariate polynomial with a
 # single variable
function isolate(f; prec = 32)
    total_degree(f) == 0 && return []
    @assert is_univariate(f) "Not univariate polynomial"

    sols = real_solutions(Ideal([change_ringvar(f)]), precision=prec, interval=true, info_level=0)
    return getindex.(sols, 1)
end

# univariate isolation of roots of a bivariate polynomial f whose
# ivar-th variable is evaluated at val
function isolate_eval(f, ivar, val; prec=64)
    fev = change_ringvar(evaluate(f, [ivar], [val]))
    # fev *= fev |> coefficients .|> denominator |> lcm
    return isolate(fev, prec=prec)
end

"""
    _expand_degenerate_intervals!(xnew, prec)

Expand intervals where lower == upper to avoid zero-width boxes.
"""
function _expand_degenerate_intervals!(xnew, prec)
    d = 1 // (ZZ(1) << prec)
    for i in eachindex(xnew)
        if xnew[i][1] == xnew[i][2]
            x = xnew[i][1]
            xnew[i] = [x - d, x + d]
        end
    end
    return xnew
end

# TODO: add the singularity check here, probably just by calling a function
# make sure we refine only the boxes that need to be refined
function  compute_crit_boxes(params, LBcrit, precx, v, max_attempts=5)
    compt = 0
    xcrit = Dict()
    while compt <= max_attempts
        try
            v > 0 && compt > 0 && println("ComputeBoxes: refine x-precision to $precx bits")
            # TODO : check that no overlap between different isolations
            xcrit = Dict(i=> reduce(vcat, [isolate(pp, prec=precx) for pp in p[1]]) for (i, p) in params)
            # Enlarge exact isolating root intervals
            for xvals in values(xcrit)
                _expand_degenerate_intervals!(xvals, precx)
            end
            arbField = ArbField(precx)
            for i in keys(xcrit)
                xvals = xcrit[i]
                ycrit = [
                    evaluate_arb(params[i][2], params[i][3], xc, arbField)
                    for xc in xvals
                ]
                _expand_degenerate_intervals!(ycrit, precx)
                @inbounds LBcrit[i] = [
                    [xvals[j], ycrit[j]] for j in eachindex(ycrit)
                ]
            end
        catch e
            precx *= 2
            v > 1 && println("Error in isolating critical boxes: ", e)
            compt += 1
        else
            break
        end
    end
    # TODO: is that a relevant upper bound?
    if compt > max_attempts
        error("Problem in isolating critical boxes")
    end

    return precx
end

function compute_singular_boxes(f, LBcrit, params, precx)
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
                    precx = compute_crit_boxes(params, LBcrit, precx, v-1)
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
    return precx
end

###################################3
# File: critical_boxes.jl

"""
    _compute_boxes_for_index!(i, params, LBcrit, precx)

Compute isolating boxes for a single index i.
"""
function _compute_boxes_for_index!(i, params, LBcrit, precx)
    xvals = reduce(vcat, (isolate(pp, prec=precx) for pp in params[i][1]))
    _expand_degenerate_intervals!(xvals, precx)

    arbField = ArbField(precx)

    yvals = [
        evaluate_arb(params[i][2], params[i][3], xc, arbField)
        for xc in xvals
    ]
    _expand_degenerate_intervals!(yvals, precx)

    @inbounds LBcrit[i] = [[xvals[j], yvals[j]] for j in eachindex(xvals)]
end


"""
    _needs_refinement(i, f, LBcrit, precx)

Check if boxes at index i require refinement.
"""
function _needs_refinement(i, f, LBcrit, precx)
    arbField = ArbField(precx)
    m = i <= 1 ? 2 : i
    # TODO: input, compute once for all i
    Lfyk = diff_list(f, 2, max(i, 2))

    for box in LBcrit[i]
        pcrit = [rat_to_arb(c, arbField) for c in box]

        if contains_zero(evaluate(Lfyk[m+1], pcrit))
            return true
        end
    end
    return false
end


"""
    compute_crit_and_singular_boxes!(f, params, LBcrit, precx; max_attempts=5, v=0)

Compute and refine critical boxes per index independently.
"""
function insulate_crit_boxes(f, params, precx; max_attempts=5, v=0)
    LBcrit, Lprecx = Dict(), Dict()

    for i in keys(params)
        attempts, precxi = 0, precx
        while attempts <= max_attempts
            try
                attempts > 0 && v > 0 &&
                    println("Refining index $i -> prec = $precxi")

                _compute_boxes_for_index!(i, params, LBcrit, precxi)

                # singularity check
                if ! needs_refinement(i, f, LBcrit, precxi, params)
                    break
                end

                # refine ONLY this index
                precxi *= 2
                attempts += 1

            catch e
                precxi *= 2
                attempts += 1
                v > 1 && println("Error at index $i: ", e)
            end
        end
        attempts > max_attempts &&
            error("Failed to isolate boxes for index $i")

        Lprecx[i] = precxi
    end

    return LBcrit, Lprecx
end
#######################################################

# Helper function to process a single edge of the box
function _isolate_box_edge(f, fixed_dim, fixed_val, target_interval, initial_prec, v; max_retries=5)
    prec = initial_prec

    for _ in 1:max_retries
        roots = isolate_eval(f, fixed_dim, fixed_val, prec=prec)
        valid_indices = Int[]
        overlap_found = false

        for (j, l) in pairs(roots)
            if in_inter(l, target_interval)
                push!(valid_indices, j)
            elseif overlap_inter(l, target_interval)
                prec *= 2
                v > 0 && println("IntersectBox: increase precision to $prec bits")
                overlap_found = true
                break # Break inner for-loop to restart with higher precision
            end
        end

        # If we successfully processed all roots without overlap, return the results
        if !overlap_found
            return (roots, valid_indices)
        end
    end

    error("Problem when isolating on one side of a box after $max_retries attempts")
end

# Function that computes the intersection of a curve with a box by isolating roots on the edges
function intersect_box(f, B; prec=100, v=0)
    # B[1] represents x-bounds, B[2] represents y-bounds

    # Evaluate horizontal edges (fix y to B[2][1] and B[2][2], target x-interval B[1])
    edge_y1 = _isolate_box_edge(f, 2, B[2][1], B[1], prec, v)::Side
    edge_y2 = _isolate_box_edge(f, 2, B[2][2], B[1], prec, v)::Side

    # Evaluate vertical edges (fix x to B[1][1] and B[1][2], target y-interval B[2])
    edge_x1 = _isolate_box_edge(f, 1, B[1][1], B[2], prec, v)::Side
    edge_x2 = _isolate_box_edge(f, 1, B[1][2], B[2], prec, v)::Side

    # Return as a 4-Tuple of 2-Tuples
    return (edge_y1, edge_y2, edge_x1, edge_x2)
end


"""
    _update_LB_first_axis!(LB, xnew)

Replace first coordinate of each box in LB with corresponding xnew interval.
"""
function _update_LB_first_axis!(LB, xnew)
    @inbounds for i in eachindex(LB)
        LB[i][1] = xnew[i]
    end
    return LB
end

"""
    refine_xboxes(f, LB, prec)

Refine LB along first axis using roots of polynomials in F
The order in F is assumed to match the one in LB
"""
function refine_xboxes(F::Vector{T} where T<:Union{PolyRingElem, MPolyRingElem}, LB, prec)
    # Concatenate all roots from isolating each polynomial in F
    xnew = reduce(vcat, [isolate(f, prec=prec) for f in F])

    _expand_degenerate_intervals!(xnew, prec)
    return _update_LB_first_axis!(LB, xnew)
end

## TODO : Refine only the intervals that need to be refined? only the factor that defines the intervals


function intersect_vertical_boxes(LBcrit, params, Lprecx)
    LPCside = Dict{Int,Any}()
    ndig = maximum([ndigits(length(LBcrit[i])) for i in eachindex(LBcrit)])
    for i in eachindex(LBcrit)
        LPCside[i] = Array{Any}(undef, length(LBcrit[i]))
        # Data printing part
        ndigi = ndigits(length(LBcrit[i]))
        Ptype = (i > length(LBcrit)-length(C)) ? "Pcon" : "mult"
        ###
        precxtmp = Lprecx[i]
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
        Lprecx[i] = precxtmp
    end
    LnPCside = Dict(i => [[length(indI) for (L, indI) in PB] for PB in LPCside[i]] for i in eachindex(LPCside))

    # Update extreme boxes
    if haskey(LBcrit, 1)
        for j in eachindex(LBcrit[1])
            # If the curve does not intersect the box only on vertical sides
            if !(LnPCside[1][j][1:2] == [0, 0])
                PCside, nPCside = LPCside[1][j], LnPCside[1][j]
                # Get the left and right sides
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
    return LPCside, LnPCside
end

###############################



function intersect_vertical_boxes_bis(LBcrit, params, Lprecx; max_attempts=5, v=0)

    LPCside = Dict{Int, Boxes}()

    for (i, boxes) in pairs(LBcrit)
        prec = Lprecx[i]
        n = length(boxes)

        sides_vec = Vector{Sides}(undef, n)

        for attempt in 1:max_attempts
            refined = false

            @inbounds for j in 1:n
                v > 0 && print("i=$i ($j/$n)\r")

                pc = intersect_box(f, boxes[j], prec=prec, v=v)::Sides

                c1 = length(pc[1][2])
                c2 = length(pc[2][2])
                c3 = length(pc[3][2])
                c4 = length(pc[4][2])

                if (i == 1 && (c1 + c2 + c3 + c4 > 2)) ||
                   (i != 1 && (c1 != 0 || c2 != 0))

                    prec *= 2
                    v > 0 && println("\nRefine -> prec = $prec")

                    refine_xboxes(params[i][1], boxes, prec)

                    refined = true
                    break
                end

                sides_vec[j] = pc
            end

            refined || break
            attempt == max_attempts &&
                error("Problem in computing intersections (i=$i)")
        end

        LPCside[i] = sides_vec
        Lprecx[i] = prec

        v > 0 && println()
    end

    # ---- Counts (type-stable) ----
    LnPCside = Dict{Int, Vector{NTuple{4, Int}}}()

    for (i, vec) in pairs(LPCside)
        LnPCside[i] = [
            (
                length(pc[1][2]),
                length(pc[2][2]),
                length(pc[3][2]),
                length(pc[4][2])
            )
            for pc in vec
        ]
    end

    # ---- Extremal case (i == 1) ----
    if haskey(LBcrit, 1)
        boxes = LBcrit[1]
        pcs   = LPCside[1]
        cnts  = LnPCside[1]

        @inbounds for j in eachindex(boxes)
            c1, c2, c3, c4 = cnts[j]

            (c1 == 0 && c2 == 0) && continue

            pc = pcs[j]
            ycrit = boxes[j][2]

            I1 = pc[3][1] # left
            I2 = pc[4][1] # right

            # Orientation of the curve left(s=1)/right(s=2)
            s = length(I1) >= length(I2) ? 1 : 2

            # ---- intersection on bottom ----
            if c1 == 1
                # Attach to the max intersection below the bow on the correct side
                yinf = maximum(i for (i, yy) in pairs((s == 1 ? I1 : I2)) if yy[1] < ycrit[1])

                push!(pc[s+2][2], yinf)
                pc[1][2] = Int[]

                cnts[j] = (
                    0,
                    c2,
                    s == 1 ? c3 + 1 : c3,
                    s == 2 ? c4 + 1 : c4
                )

                c1, c2, c3, c4 = cnts[j] # for next step
            end

            # ---- intersection on top ----
            if c2 == 1
                ymax = minimum(i for (i, yy) in pairs((s == 1 ? I1 : I2)) if yy[2] > ycrit[2])

                push!(pc[s+2][2], ymax)
                pc[2][2] = Int[]

                cnts[j] = (
                    c1,
                    0,
                    s == 1 ? c3 + 1 : c3,
                    s == 2 ? c4 + 1 : c4
                )
            end
        end
    end

    return LPCside, LnPCside
end