# for convenience
function affine_cell(F::Vector{T},
                     H::Vector{T}) where {T <: MPolyRingElem}
    return LocClosedSet{T}(F,H)
end

function num_eqns(X::LocClosedSet)
    length(findall((!).(X.eqns_is_red)))
end

# for displaying locally closed sets
function Base.show(io::IO, lc::LocClosedSet)
    string_rep = "V("
    for (i, f) in enumerate(lc.eqns)
        if i != length(lc.eqns)
            string_rep *= "$f, "
        else
            string_rep *= "$(f)) \\ "
        end
    end
    string_rep *= "V("
    for (i, f) in enumerate(lc.ineqns)
        if i != length(lc.ineqns)
            string_rep *= "($f)*"
        else
            string_rep *= "($f)"
        end
    end
    string_rep *= ")"
    print(io, string_rep)
end

function ring(X::LocClosedSet)
    return parent(first(X.eqns))
end

function codim(X::LocClosedSet)
    mis = first(max_ind_sets(X.gb))
    return length(findall(b -> !b, mis))
end

function is_empty_set(X::LocClosedSet)
    R = ring(X)
    if one(R) in X.gb
        return true
    else
        gb2 = saturate(X.gb, last(gens(R)))
        return one(R) in gb2
    end
    return false
end

function add_equation!(X::LocClosedSet, f::MPolyRingElem)
    @info "adding equation"
    push!(X.eqns, f)
    push!(X.eqns_is_red, false)
    X.gb = saturate(vcat(X.gb, [f]), X.ineqns)
end

function add_inequation!(X::LocClosedSet, h::MPolyRingElem; method=:sat)
    # @info "adding inequation"
    push!(X.ineqns, h)
    X.gb = method == :sat ? saturate(X.gb, h) : quotient(X.gb, h)
end

function add_inequation(X::LocClosedSet, h::MPolyRingElem; kwargs...)
    if isone(h)
        return X
    end
    Y = deepcopy(X)
    add_inequation!(Y, h; kwargs...)
    return Y
end

function add_inequations(X::LocClosedSet{P}, H::Vector{P}) where P
    Htil = filter(!isone, H)
    @info "adding ineqns $((Nemo.leading_monomial).(H))"
    Y = deepcopy(X)
    append!(Y.ineqns, Htil)
    Y.gb = saturate(Y.gb, Htil)
    return Y
end

function hull(X::LocClosedSet, g::MPolyRingElem; method = :sat)
    @info "hull for lm $(Nemo.leading_monomial(g))"
    gb = X.gb
    col_gb = method == :sat ? saturate(gb, g) : quotient(gb, g)
    if one(ring(X)) in col_gb
        return [X]
    end
    sort!(col_gb, by = p -> total_degree(p))
    H_rand = filter(!iszero, my_normal_form(col_gb, gb))
    if isempty(H_rand)
        @info "regular intersection in hull"
        return typeof(X)[]
    end
    res = remove(X, H_rand)
    # res = remove_with_tree(X, H)
    return res
end

# function remove(X::LocClosedSet,
#                 H::Vector{<:MPolyRingElem})
    
#     res = typeof(X)[]
#     isempty(H) && return res
#     h = first(H)
#     cells = remove(X, H[2:end])
#     if iszero(my_normal_form([h], X.gb))
#         return cells
#     else
#         @info "inequation for lm $(Nemo.leading_monomial(h))"
#         Y = add_inequation(X, h)
#         if is_empty_set(Y)
#             @info "empty set in remove"
#             return cells
#         end
#         push!(res, Y)
#         for Z in cells
#             cells2 = hull(Z, h)
#             append!(res, cells2)
#         end
#         return res
#     end
#     return res
# end

function remove(X::LocClosedSet,
                H::Vector{<:MPolyRingElem},
                ineq::MPolyRingElem=one(ring(X)))

    res = typeof(X)[]
    isempty(H) && return res

    h = first(H)
    @info "adding ineqn $((Nemo.leading_monomial(h), Nemo.leading_monomial(ineq)))"
    Y = add_inequation(X, h*ineq)
    if is_empty_set(Y)
        return remove(X, H[2:end], ineq)
    end
    push!(res, Y)
    G = filter(!iszero, my_normal_form(random_lin_combs(Y.gb) .* ineq, X.gb))
    for htil in H[2:end]
        @info "taking hulls for $(Nemo.leading_monomial(h))"
        new_comps = remove(X, G, htil)
        append!(res, new_comps)
    end
    return res
end

# --- helper functions --- #

function saturate(F::Vector{P}, nz::P) where {P <: MPolyRingElem}
    return saturate(F, [nz])
end

function saturate(F::Vector{P}, nzs::Vector{P}) where {P <: MPolyRingElem}
    R = parent(first(F))
    S, vars = polynomial_ring(base_ring(R), vcat(["t$i" for i in 1:length(nzs)], ["x$i" for i in 1:nvars(R)]),
                              ordering = :degrevlex)
    Fconv = [convert_poly_to_t_ring(f, S) for f in F]

    for (i, h) in enumerate(nzs)
        ti = vars[i]
        push!(Fconv, ti*convert_poly_to_t_ring(h, S)-1)
    end

    gb = groebner_basis(Ideal(Fconv), complete_reduction = true,
                        eliminate = length(nzs))

    # convert back to original ring
    res = Vector{P}(undef, length(gb))

    for (i, p) in enumerate(gb)
        res[i] = convert_to_orig_ring(p, R) 
    end
    return res
end

function quotient(F::Vector{P}, nz::P) where {P <: MPolyRingElem}
    return quotient(F, [nz])
end

function quotient(F::Vector{P}, nzs::Vector{P}) where {P <: MPolyRingElem}
    R = parent(first(F))
    S, vars = polynomial_ring(base_ring(R), vcat(["t$i" for i in 1:length(nzs)], ["x$i" for i in 1:nvars(R)]),
                              ordering = :degrevlex)
    Fconv = [convert_poly_to_t_ring(f, S) for f in F]

    for (i, h) in enumerate(nzs)
        ti = vars[i]
        Fconv .*= ti
        push!(Fconv, (ti-1)*convert_poly_to_t_ring(h, S))
    end

    gb = groebner_basis(Ideal(Fconv), complete_reduction = true,
                        eliminate = length(nzs))

    # convert back to original ring
    res = Vector{P}(undef, length(gb))

    for (i, p) in enumerate(gb)
        res[i] = divides(convert_to_orig_ring(p, R), prod(nzs))[2]
    end
    return res
end

# assumes H is sorted by degree
function random_lin_combs(H::Vector{P}) where {P <: MPolyRingElem}
    res = P[]
    chr = characteristic(base_ring(first(H)))
    curr_deg = total_degree(first(H))
    curr_pol = zero(parent(first(H)))
    for h in H
        if total_degree(h) == curr_deg
            curr_pol += rand(1:chr-1)*h
        else
            push!(res, curr_pol)
            curr_pol = h
            curr_deg = total_degree(h)
        end
    end
    push!(res, curr_pol)
    return res
end

function random_lin_comb(F::Vector{P}) where {P <: MPolyRingElem}
    R = parent(first(F))
    res = zero(R)
    chr = characteristic(R)
    for f in F
        res += rand(1:chr-1)*f
    end
    return res
end

function convert_poly_to_t_ring(f::P, S::MPolyRing) where {P <: MPolyRingElem}
    ctx = MPolyBuildCtx(S)
    R = parent(f)
    nts = nvars(S) - nvars(R)
    for (e, c) in zip(exponent_vectors(f), coefficients(f))
        enew = vcat(zeros(Int, nts), e)
        push_term!(ctx, c, enew)
    end
    return finish(ctx)
end

function convert_to_orig_ring(f::P, R::MPolyRing) where {P <: MPolyRingElem}
    ctx = MPolyBuildCtx(R)
    S = parent(f)
    nts = nvars(S) - nvars(R)
    for (e, c) in zip(exponent_vectors(f), coefficients(f))
        push_term!(ctx, c, e[nts+1:end])
    end
    return finish(ctx)
end

function max_ind_sets(gb::Vector{P}) where {P <: MPolyRingElem}
    R = parent(first(gb))
    res = [trues(ngens(R))]

    lms = (Nemo.leading_monomial).(gb)
    for lm in lms
        to_del = Int[]
        new_miss = BitVector[]
        for (i, mis) in enumerate(res)
            nz_exps_inds = findall(e -> !iszero(e),
                                   first(Nemo.exponent_vectors(lm)))
            ind_var_inds = findall(mis)
            if issubset(nz_exps_inds, ind_var_inds)
                for j in nz_exps_inds
                    new_mis = copy(mis)
                    new_mis[j] = false
                    push!(new_miss, new_mis)
                end
                push!(to_del, i)
            end
        end
        deleteat!(res, to_del)
        append!(res, new_miss)
        unique!(res)
    end

    max_length = maximum(mis -> length(findall(mis)), res)
    filter!(mis -> length(mis) != max_length, res)
    return res
end
    

# for debugging
function check_decomp(F::Vector{P}, Xs::Vector{<:LocClosedSet}) where {P <: MPolyRingElem}
    println("checking decomp")
    gb_ch = F
    for X in Xs
        g = random_lin_comb(X.gb)
        gb_ch = saturate(gb_ch, g)
    end
    R = parent(first(F))
    gb_ch = saturate(gb_ch, last(gens(R)))
    return one(parent(first(F))) in gb_ch
end

# EXPERIMENTAL REMOVE

# mutable struct NZNode{T}
#     p::T
#     parent::Union{Nothing, NZNode{T}}
#     children::Vector{NZNode{T}}
# end

# AbstractTrees.children(nd::NZNode) = nd.children
# AbstractTrees.nodevalue(nd::NZNode) = Nemo.leading_monomial(nd.p)
# AbstractTrees.ParentLinks(::Type{<:NZNode}) = StoredParents()
# AbstractTrees.parent(nd::NZNode) = nd.parent
# AbstractTrees.NodeType(::Type{<:NZNode}) = HasNodeType()
# AbstractTrees.nodetype(::Type{<:NZNode{T}}) where T = NZNode{T}

# function add_child!(nd::NZNode{T}, p::T) where T
#     ch = NZNode(p, nd, NZNode{T}[])
#     push!(nd.children, ch)
# end

# function right_siblings(nd::NZNode{T}) where T
#     pr = AbstractTrees.parent(nd)
#     if isnothing(pr)
#         return NZNode{T}[]
#     else
#         ch_idx = findfirst(ch -> ch == nd, children(pr))
#         return children(pr)[ch_idx+1:end]
#     end
# end

# function right_sibling(nd::NZNode)
#     pr = AbstractTrees.parent(nd)
#     if isnothing(pr)
#         return nothing
#     else
#         ch_idx = findfirst(ch -> ch == nd, children(pr))
#         if ch_idx == length(children(pr))
#             return nothing
#         else
#             return children(pr)[ch_idx+1]
#         end
#     end
# end

# function get_pols(nd::NZNode)
#     res = [nd.p]
#     pr = AbstractTrees.parent(nd) 
#     while !isnothing(pr)
#         p = pr.p
#         !isone(p) && push!(res, pr.p)
#         pr = AbstractTrees.parent(pr)
#     end
#     return res
# end

# function next_leaf(nd::NZNode)
#     cr_node = nd
#     while !isnothing(cr_node) && !isroot(cr_node)
#         n = right_sibling(cr_node)
#         if isnothing(n)
#             cr_node = AbstractTrees.parent(cr_node)
#         else
#             # descendleft is not allowed to return nothing
#             return descendleft(n)
#         end
#     end
#     nothing
# end

# function remove_with_tree(X::LocClosedSet{P}, H::Vector{P}) where P
#     res = typeof(X)[]

#     R = ring(X)
#     root = NZNode(one(R), nothing, NZNode{P}[])

#     H_rand = random_lin_combs(H)
#     for h in H_rand
#         add_child!(root, h)
#     end

#     cr_node = first(children(root))
#     while !isnothing(cr_node)
#         pls = get_pols(cr_node)
#         Y = add_inequations(X, pls)
#         if !is_empty_set(Y)
#             push!(res, Y)

#             new_ineqns = random_lin_combs(filter(!iszero, my_normal_form(Y.gb, X.gb)))
#             sort!(new_ineqns, by = p -> total_degree(p))
#             for sbl in right_siblings(cr_node)
#                 if !isempty(children(sbl))
#                     error(":(")
#                 end
#                 for g in new_ineqns
#                     add_child!(sbl, g)
#                 end
#             end
#         end
#         cr_node = next_leaf(cr_node)
#     end

#     return res
# end
