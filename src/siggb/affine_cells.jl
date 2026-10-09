function loc_closed_set(seq::Vector{T}) where {T<:MPolyRingElem}
    @assert !isempty(seq) "cannot construct affine cell from no equations."
    R = parent(first(seq))
    codim_upper_bound = min(length(seq), ngens(R) - 1)
    gb = saturate(seq, last(gens(R)))
    # index 0 indicates that no rational reconstruction will need to be performed
    return LocClosedSet([(f, 0) for f in seq], codim_upper_bound, [gb], [Int[]])
end

# basic data
function num_eqns(X::LocClosedSet)
    return length(X.seq)
end

function ring(X::LocClosedSet)
    return parent(first(first(X.seq)))
end

# queries
function is_empty_set(X::LocClosedSet)
    R = ring(X)
    if isempty(X.gbs)
        return true
    elseif all(gb -> one(R) in gb, X.gbs)
        return true
    elseif all(isempty, X.gbs)
        return true
    end
    return false
end

function is_equidimensional(X::LocClosedSet)
    return all(gb -> codim(gb) == X.codim_upper_bound, X.gbs)
end

# core functions
function add_to_output!(res::Vector{C},
                        X::C) where C<:LocClosedSet

    if is_empty_set(X)
        @info "empty component"
        return true
    elseif is_equidimensional(X)
        @info "component with $(num_eqns(X)) eqns done, codimension $(X.codim_upper_bound)"
        push!(res, X)
        return true
    end
    return false
end

function add_inequation!(X::LocClosedSet, h::P, r::Registry;
                         known_zds=P[]::Vector{P}) where P
    
    h *= leading_coefficient(h)^(-1)
    ri = update_registry!(r, h)
    for (i, gb) in enumerate(X.gbs)
        X.gbs[i] = saturate(vcat(gb, known_zds), h)
        push!(X.ineqns[i], ri)
    end
    return ri
end

function add_inequation(X::LocClosedSet, h::MPolyRingElem, r::Registry)
    Y = deepcopy(X)
    if isone(h)
        return Y, 0
    end
    ri = add_inequation!(Y, h, r)
    return Y, ri
end

# which equations are hull equations needs to be managed outside of this function
function split(X::LocClosedSet, g::MPolyRingElem, r::Registry)
    tim = @elapsed X_min_g, ri = add_inequation(X, g, r)
    @info "initial saturation time $(tim)"

    X_hull_g = deepcopy(X)

    col_gbs = X_min_g.gbs
    hull_gbs = Vector{typeof(g)}[]
    ineqns = Vector{Int}[]
    R = ring(X)
    for (X_gb, X_gb_ineqns, col_gb) in zip(X.gbs, X.ineqns, col_gbs)
        if one(R) in col_gb
            @info "equation vanishes on one GB"
            tim = @elapsed new_gb = saturate(vcat(X_gb, [g]), last(gens(R)))
            @info "adding equation time $(tim)"
            push!(hull_gbs, new_gb)
            push!(ineqns, X_gb_ineqns)
            continue
        end
        sort(col_gb, by = p -> total_degree(p))
        H_rand = filter(!iszero, normal_form(random_lin_combs(col_gb, r.rand_coeffs), X_gb))
        isempty(H_rand) && continue
        gbsineqns = remove!(X_gb, H_rand, r, known_eqns = [g])
        for (gb, gb_ineqns) in gbsineqns
            push!(hull_gbs, gb)
            push!(ineqns, vcat(X_gb_ineqns, gb_ineqns))
        end
    end
    X_hull_g.gbs = hull_gbs
    X_hull_g.ineqns = ineqns

    todel = findall(gb -> one(R) in gb, X_min_g.gbs)
    deleteat!(X_min_g.gbs, todel)
    deleteat!(X_min_g.ineqns, todel)

    todel = findall(gb -> one(R) in gb, X_hull_g.gbs)
    deleteat!(X_hull_g.gbs, todel)
    deleteat!(X_hull_g.ineqns, todel)

    return X_hull_g, X_min_g, ri
end

function remove!(gb::Vector{P},
                 H::Vector{P},
                 r::Registry;
                 known_eqns::Vector{P}=P[]) where P

    res = Tuple{Vector{P}, Vector{Int}}[]
    isempty(H) && return res

    R = parent(first(gb))
    h = first(H)
    tim = @elapsed gb1 = saturate(vcat(gb, known_eqns), h)
    @info "remove time $(tim) for degree $(total_degree(h))"
    if one(R) in gb1
        @info "is empty"
        return remove!(gb, H[2:end], r, known_eqns=known_eqns)
    end
    h *= leading_coefficient(h)^(-1)
    ri = update_registry!(r, h)
    push!(res, (gb1, [ri]))
    tim1 = @elapsed G = filter(!iszero,
                               normal_form(random_lin_combs(gb1, r.rand_coeffs), gb))
    if isempty(G)
        return res
    end
    g_rand = random_lin_comb(G, r.rand_coeffs)
    rem_rest = H[2:end]
    tim2 = @elapsed filter!(h -> !iszero(normal_form(h*g_rand, gb)), rem_rest)
    @info "normal forms computed in $(tim1 + tim2)"
    gbsineqns2 = remove!(gb, rem_rest, r, known_eqns = known_eqns)
    for (gb2, ineqns2) in gbsineqns2
        gbs3ineqns = remove!(gb2, G, r, known_eqns = vcat(known_eqns, [h]))
        append!(res, [(gb3, vcat(ineqns2, ineqns3)) for (gb3, ineqns3) in gbs3ineqns])
    end
    return res
end

# -------------------------------------- #
# --- user level locally closed sets --- #
# -------------------------------------- #

function Base.show(io::IO, ::MIME"text/plain", X::LocallyClosedSet)
    str = "V("
    for (i, eqn) in enumerate(X.eqns)
        str *= i == length(X.eqns) ? "$(eqn))" : "$(eqn), "
    end

    if !isempty(X.ineqns)
        str *= " \\ V("
        for (i, ineqn) in enumerate(X.ineqns)
            ineqn_repr = if isone(length(ineqn)) || isone(length(X.ineqns))
                "$(ineqn)"
            else
                "($(ineqn))"
            end
            str *= i == length(X.ineqns) ? "$(ineqn_repr))" : "$(ineqn_repr) * "
        end
    end

    print(io, str)
end

@doc Markdown.doc"""
    equations(X::LocallyClosedSet)

Given a locally closed set `X` of the form $V(F) \ V(g_1 \cdot \dots \cdot g_r)$,
return the list of polynomials $F$.
"""
equations(X::LocallyClosedSet) = X.eqns

@doc Markdown.doc"""
    inequations(X::LocallyClosedSet)

Given a locally closed set `X` of the form $V(F) \ V(g_1 \cdot
\dots \cdot g_r)$, return the list of polynomials $g_1, \dots g_r$.
"""
inequations(X::LocallyClosedSet) = X.ineqns

@doc Markdown.doc"""
    dimension(X::LocallyClosedSet)

Return the dimension of `X`.
"""
dimension(X::LocallyClosedSet) = X.dim

@doc Markdown.doc"""
    Ideal(X::LocallyClosedSet)

Return a polynomial ideal `I` whose zero locus coincides with the
Zariski closure of `X`.

**Note**: If no ideal `I` is known, this will require a Gröbner basis
computation. The output of this function is cached.
"""
function Ideal(X::LocallyClosedSet)
    !ismissing(X.ideal) && return X.ideal
    gb = saturate(X.eqns, X.ineqns)
    I = Ideal(gb)
    I.gb[0] = gb
    return I
end

# convert internal cells to user level output cells
function get_output_cells(cell::LocClosedSet,
                          R::FqMPolyRing,
                          r::ModularRegistry)

    res = LocallyClosedSet{FqMPolyRingElem}[]
    eqns = _dehomogenize(first.(cell.seq), R)
    dim = ngens(R) - cell.codim_upper_bound
    for (gb, ineqninds) in zip(cell.gbs, cell.ineqns)
        gb_dehom = _dehomogenize(gb, R)
        ls = if isempty(ineqninds)
            LocallyClosedSet(eqns, dim)
        else
            ineqns = unique(_dehomogenize(get_pols(r, ineqninds), R))
            LocallyClosedSet(eqns, ineqns, dim)
        end
        I = Ideal(gb_dehom)
        I.gb[0] = gb_dehom
        ls.ideal = I
        push!(res, ls)
    end
    return res
end

# convert internal cells to user level output cells
function get_output_cells(cell::LocClosedSet,
                          R::QQMPolyRing,
                          input_eqns::Vector{QQMPolyRingElem},
                          r::ReconstructRegistry)

    res = LocallyClosedSet{QQMPolyRingElem}[]
    S = ring(cell)
    eqns = QQMPolyRingElem[]
    for (f, i) in cell.seq
        if iszero(i)
            inp_ind = findfirst(g -> reduce_mod_p(g, S) == f, input_eqns)
            push!(eqns, input_eqns[inp_ind])
        else
            push!(eqns, get_pol(r, i))
        end
    end
    eqns = _dehomogenize(eqns, R)
    dim = ngens(R) - cell.codim_upper_bound
    for ineqninds in cell.ineqns
        if isempty(ineqninds)
            push!(res, LocallyClosedSet(eqns, dim))
        else
            ineqns = unique(_dehomogenize(get_pols(r, ineqninds), R))
            push!(res, LocallyClosedSet(eqns, ineqns, dim))
        end
    end
    return res
end

# ------------------------ #
# --- helper functions --- #
# ------------------------ #

function saturate(F::Vector{P}, nz::P) where {P <: MPolyRingElem}
    return saturate(F, [nz])
end

function saturate(F::Vector{P}, nzs::Vector{P}) where {P <: MPolyRingElem}
    R = parent(first(F))
    S, vars = polynomial_ring(base_ring(R), vcat(["t$i" for i in 1:length(nzs)], ["x$i" for i in 1:nvars(R)]),
                              internal_ordering = :degrevlex)
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
                              internal_ordering = :degrevlex)
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

# kept small: the reconstructed polynomials inherit the height of these
# coefficients, so larger ones would cost extra primes
const RAND_COEFF_BOUND = 10000

reset_rand_coeffs!(rc::RandCoeffs) = rc.ind = 1

function next_rand_coeff!(rc::RandCoeffs, K)
    i = rc.ind
    while length(rc.coeffs) < i
        c = rand(1:RAND_COEFF_BOUND)
        push!(rc.coeffs, rand(Bool) ? c : -c)
    end
    rc.ind = i + 1
    return K(rc.coeffs[i])
end

# assumes H is sorted by degree
function random_lin_combs(H::Vector{P}, rc::RandCoeffs) where {P <: MPolyRingElem}
    res = P[]
    K = base_ring(first(H))
    curr_deg = total_degree(first(H))
    curr_pol = zero(parent(first(H)))
    for h in H
        if total_degree(h) == curr_deg
            curr_pol += next_rand_coeff!(rc, K)*h
        else
            push!(res, curr_pol)
            curr_pol = h
            curr_deg = total_degree(h)
        end
    end
    push!(res, curr_pol)
    return res
end

function random_lin_comb(F::Vector{P}, rc::RandCoeffs) where {P <: MPolyRingElem}
    R = parent(first(F))
    res = zero(R)
    K = base_ring(R)
    for f in F
        res += next_rand_coeff!(rc, K)*f
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
    filter!(mis -> length(findall(mis)) == max_length, res)
    return res
end

function codim(gb::Vector{P}) where P
    miss = max_ind_sets(gb)
    cd = length(findall(b -> !b, first(miss)))
    return cd
end

function variety_string_rep(F::Vector{<:MPolyRingElem};
                            sep = ", ", lpar = "", rpar = "")
    string_rep = "V("
    for (i, f) in enumerate(F)
        if isone(i)
            string_rep *= (lpar * "$(f)" * rpar)
        else
            string_rep *= (sep * lpar * "$(f)" * rpar)
        end
    end
    string_rep *= ")"
    return string_rep
end

# check if idls provides ideal-theoretic decomposition of radical of I
function _check_decomp(I::IDL, Xs::Vector{<:LocallyClosedSet}) where {IDL <: Ideal}
    gb_ch = I.gens
    rc = RandCoeffs()
    for X in Xs
        idl = Ideal(X)
        g = random_lin_comb(idl.gens, rc)
        gb_ch = saturate(gb_ch, g)
    end
    R = parent(I)
    return one(R) in gb_ch
end
