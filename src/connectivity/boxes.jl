in_inter(I, J) =  J[1] <= I[1] && I[2] <= J[2]
overlap_inter(I,J) = max(I[1], J[1]) <= min(I[2], J[2])

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
                v > 0 && println("Increase precision to ", prec)
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

function refine_xboxes(f::T where T<:Union{PolyRingElem, MPolyRingElem}, LB, prec)
    # Refine LB along first axis, being roots of f
    xnew = isolate(f, prec=prec)
    for i in eachindex(xnew)
        if xnew[i][1]==xnew[i][2]
            xnew[i] = [xnew[i][1]-1//ZZ(1)<<prec, xnew[i][1]+1//ZZ(1)<<prec]
        end
    end
    for i in eachindex(LB)
		  LB[i] = [ xnew[i], LB[i][2] ]
    end
end

function refine_xboxes(F::Vector{T} where T<:Union{PolyRingElem, MPolyRingElem}, LB, prec)
    # Refine LB along first axis, being roots of the polynomial in F
    # The order in F is assumed to match the one in LB
    xnew = reduce(vcat, [isolate(f, prec=prec) for f in F])
    for i in eachindex(xnew)
        if xnew[i][1]==xnew[i][2]
            xnew[i] = [xnew[i][1]-1//ZZ(1)<<prec, xnew[i][1]+1//ZZ(1)<<prec]
        end
    end
    for i in eachindex(LB)
		  LB[i] = [ xnew[i], LB[i][2] ]
    end
end

function refine_ith_xboxes(f, LB, i, prec)
    # Refine only LB[i] along first axis, being i-th root of f
    xnew = isolate(f, prec=prec)
	LB[i] = [ xnew[i], LB[i][2] ]
end
