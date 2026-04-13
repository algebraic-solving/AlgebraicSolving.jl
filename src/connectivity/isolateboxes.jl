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
    edge_y1 = _isolate_box_edge(f, 2, B[2][1], B[1], prec, v)
    edge_y2 = _isolate_box_edge(f, 2, B[2][2], B[1], prec, v)

    # Evaluate vertical edges (fix x to B[1][1] and B[1][2], target y-interval B[2])
    edge_x1 = _isolate_box_edge(f, 1, B[1][1], B[2], prec, v)
    edge_x2 = _isolate_box_edge(f, 1, B[1][2], B[2], prec, v)

    # Return as a typed array of Tuples
    return [edge_y1, edge_y2, edge_x1, edge_x2]
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

"""
    refine_ith_xboxes(f, LB, i, prec)

Refine only the i-th box along first axis using i-th root of f.
"""
function refine_ith_xboxes(f, LB, i, prec)
    xnew = isolate(f, prec=prec)
    @inbounds LB[i] = [xnew[i], LB[i][2]]
    return LB
end