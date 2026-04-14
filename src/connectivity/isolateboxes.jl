in_inter(I, J) =  J[1] <= I[1] && I[2] <= J[2]
overlap_inter(I, J) = max(I[1], J[1]) <= min(I[2], J[2])

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
    fev = evaluate(f, [ivar], [val])
    # fev *= fev |> coefficients .|> denominator |> lcm
    fev = fev / content(fev)
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
function _needs_refinement(i, Lfyk, LBcrit, precx)
    arbField = ArbField(precx)
    m = i <= 1 ? 2 : i
    for box in LBcrit[i]
        pcrit = [rat_to_arb(c, arbField) for c in box]

        if contains_zero(evaluate(Lfyk[m+1], pcrit))
            return true
        end
    end
    return false
end

# Compute the list of the n-th first derivative of p w.r.t v
function diff_list(p, v, n)
    L = Vector{typeof(p)}(undef, n + 1)
    L[1] = p
    @inbounds for j in 2:n+1
        L[j] = derivative(L[j-1], v)
    end
    return L
end

"""
    compute_crit_and_singular_boxes(f, params, precx; max_attempts=5, v=0)

Compute insulating boxes for critical points, refining precision if necessary.
Insulating means that each box contains a single critical point and
all branches of the curve intersecting the box are incident to this
critical point.
"""
function insulate_crit_boxes(f, params, precx; max_attempts=5, v=0)
    LBcrit, Lprecx = Dict(), Dict()
    Lfyk = diff_list(f, 2, max(maximum(keys(params)), 2))

    for i in keys(params)
        attempts, precxi = 0, precx
        while attempts <= max_attempts
            try
                attempts > 0 && v > 0 &&
                    println("Refining index $i -> prec = $precxi")

                _compute_boxes_for_index!(i, params, LBcrit, precxi)

                # singularity check
                if !_needs_refinement(i, Lfyk, LBcrit, precxi)
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

# Helper function to process a single edge of the box
function _isolate_box_edge(f, fixed_dim, fixed_val, target_interval, initial_prec, v, outf; max_retries=5)
    prec = initial_prec

    for _ in 1:max_retries
        roots = isolate_eval(f, fixed_dim, fixed_val, prec=prec)
        indices_inside = Int[]
        overlap_found = false

        for (j, l) in pairs(roots)
            if in_inter(l, target_interval)
                push!(indices_inside, j)
            elseif overlap_inter(l, target_interval)
                prec *= 2
                v > 0 && println("IntersectBox: increase precision to $prec bits")
                overlap_found = true
                break
            end
        end
        # All roots without overlap on this edge
        if !overlap_found
            T = QQFieldElem
            return BoxEdge{T}(roots, indices_inside)
        end
    end
    error("Problem when isolating on one side of a box after $max_retries attempts")
end

# Function that computes the intersection of a curve with a box by isolating roots on the edges
function intersect_box(f, B, outf; prec=100, v=0)
    # Evaluate horizontal edges (fix y to B[2][1] and B[2][2], target x-interval B[1])
    edge_y1 = _isolate_box_edge(f, 2, B[2][1], B[1], prec, v, outf)
    edge_y2 = _isolate_box_edge(f, 2, B[2][2], B[1], prec, v, outf)
    # Evaluate vertical edges (fix x to B[1][1] and B[1][2], target y-interval B[2])
    edge_x1 = _isolate_box_edge(f, 1, B[1][1], B[2], prec, v, outf)
    edge_x2 = _isolate_box_edge(f, 1, B[1][2], B[2], prec, v, outf)

    # Return clean struct
    T = QQFieldElem
    return BoxIntersections{T}(edge_y1, edge_y2, edge_x1, edge_x2)
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

# Function computing the intersection of the vertical sides of the critical boxes with the curve defined by f.
# If an intersection on top/bottom side is detected, refine the width until this does not happen or max_attemps is reached.
function intersect_vertical_boxes(f, params, LBcrit, Lprecx, outf; max_attempts=5, v=0)
    # Automatically infer the point type based on outf
    T = QQFieldElem
    LPCside = Dict{Int, Vector{BoxIntersections{T}}}()
    LnPCside = Dict{Int, Vector{NTuple{4, Int}}}()

    for (i, boxes) in pairs(LBcrit)
        prec = Lprecx[i]
        n = length(boxes)
        sides_vec = Vector{BoxIntersections{T}}(undef, n)
        nsides_vec = Vector{NTuple{4, Int}}(undef,n)

        for attempt in 1:max_attempts
            refined = false

            @inbounds for j in 1:n
                v >= 0 && print("i=$i ($j/$n)\r")
                pc = intersect_box(f, boxes[j], outf, prec=prec, v=v)

                c1 = length(pc.bottom.indices_inside)
                c2 = length(pc.top.indices_inside)
                c3 = length(pc.left.indices_inside)
                c4 = length(pc.right.indices_inside)

                if (i == 1 && (c1 + c2 + c3 + c4 > 2)) ||
                   (i != 1 && (c1 != 0 || c2 != 0))
                    prec *= 2
                    v > 0 && println("\nRefine -> prec = $prec")
                    refine_xboxes(params[i][1], boxes, prec)
                    refined = true
                    break
                end

                sides_vec[j] = pc
                nsides_vec[j] = (c1,c2,c3,c4)
            end

            refined || break
            attempt == max_attempts &&
                error("Problem in computing intersections (i=$i)")
        end

        LPCside[i] = sides_vec
        LnPCside[i] = nsides_vec
        Lprecx[i] = prec
        v >= 0 && println()
    end

    # ---- Extremal case (i == 1) we vertically enlarge boxes manually -----
    # TODO: should we also modify LBcrit?
    if haskey(LBcrit, 1)
        boxes, pcs, npcs = LBcrit[1], LPCside[1], LnPCside[1]

        @inbounds for j in eachindex(boxes)
            c1, c2, c3, c4 = npcs[j]

            (c1 == 0 && c2 == 0) && continue # Nothing to do

            pc = pcs[j]
            ycrit = boxes[j][2]

            # Orientation of the curve left(s=1)/right(s=2)
            s = length(pc.left.points) >= length(pc.right.points) ? 1 : 2
            side_pts = s == 1 ? pc.left.points : pc.right.points
            side_indices = s == 1 ? pc.left.indices_inside : pc.right.indices_inside

            # ---- intersection on bottom ----
            if c1 == 1
                # Attach to the max intersection below the box on the correct side
                ind_yinf = maximum(i for (i, yy) in pairs(side_pts) if yy < ycrit[1])
                push!(side_indices, ind_yinf)
                empty!(pc.bottom.indices_inside)
                npcs[j] = (0, c2, s == 1 ? c3 + 1 : c3, s == 2 ? c4 + 1 : c4)
            end

            c1, c2, c3, c4 = npcs[j] # In case this has changed
            # ---- intersection on top ----
            if c2 == 1
                # Attach to the min intersection above the bow on the correct side
                ind_ymax = minimum(i for (i, yy) in pairs(side_pts) if yy > ycrit[2])
                push!(side_indices, ind_ymax)
                empty!(pc.top.indices_inside)
                npcs[j] = (c1, 0, s == 1 ? c3 + 1 : c3, s == 2 ? c4 + 1 : c4)
            end
        end
    end

    return LPCside, LnPCside
end