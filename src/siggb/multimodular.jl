function reduce_mod_p(f::QQMPolyRingElem, R::FqMPolyRing)

    F = base_ring(R)
    ctx = MPolyBuildCtx(R)
    for (cf, exp) in zip(coefficients(f), exponent_vectors(f))
        cfP = F(numerator(cf)) * F(denominator(cf))^(-1)
        push_term!(ctx, cfP, exp)
    end
    return finish(ctx)
end

function is_finished(r::ReconstructRegistry)
    !isempty(r.pols) && iszero(r.n_unstable)
end

function get_pol(r::ReconstructRegistry, i::Int)
    rp = r.pols[i]
    ctx = MPolyBuildCtx(r.R)
    for (cf, exp) in zip(rp.coeff_cands, rp.exps)
        push_term!(ctx, cf, Int64.(exp))
    end
    return finish(ctx)
end

function get_pols(r::ReconstructRegistry, inds::Vector{Int})
    return [get_pol(r, i) for i in inds]
end

function get_pols(r::ModularRegistry, inds::Vector{Int})
    return r.pols[inds]
end

function new_prime!(r::ReconstructRegistry, p::Integer)
    if !iszero(r.current_prime)
        push!(r.primes, r.current_prime)
        r.pprod *= r.current_prime
    end
    r.current_prime = p
    r.curr_ind = 1
end

# TODO:adjust
function does_not_match(pr::ReconstructPol, p::FqMPolyRingElem)
    ev = collect(exponent_vectors(p))
    return length(ev) != length(pr.exps) || first(ev) != first(pr.exps)
end

attempt_reconstruction(nprimes::Int) = nprimes >= 2 && ispow2(nprimes)

# Check a candidate reconstructed from pprod against a prime not dividing pprod.
function verifies_mod_p(q::QQFieldElem, r, F)
    return F(numerator(q)) == r * F(denominator(q))
end

function ReconstructPol(p::FqMPolyRingElem)
    p *= leading_coefficient(p)^(-1)
    exps = collect(exponent_vectors(p))
    mod_coeffs = (c -> lift(ZZ, c)).(collect(coefficients(p)))
    coeff_cands = (c -> QQ(c)).(mod_coeffs)
    return ReconstructPol(exps, coeff_cands, mod_coeffs, [false for _ in coeff_cands])
end

# update polynomial at currend index of registry
function update_registry!(reg::ReconstructRegistry,
                          new_pol::FqMPolyRingElem)


    ri = reg.curr_ind

    @info "updating registry at index $ri"

    if length(reg.pols) < ri
        new_rp = ReconstructPol(new_pol)
        push!(reg.pols, new_rp)
        reg.n_unstable += length(new_rp.is_stable)
        reg.curr_ind += 1
        return ri
    end

    pr = reg.pols[ri]
    does_not_match(pr, new_pol) && error("Bad prime during multi-modular computation.")

    if all(pr.is_stable)
        reg.curr_ind += 1
        return ri
    end

    pprod = reg.pprod
    curr_p = reg.current_prime
    do_reconstruct = attempt_reconstruction(length(reg.primes))
    F = base_ring(parent(new_pol))

    i = 1
    for (ccurr, cnew_fq) in zip(pr.mod_coeffs, coefficients(new_pol))
        if pr.is_stable[i]
            i += 1
            continue
        end
        if do_reconstruct
            success, new_qq_coeff = unsafe_reconstruct(ccurr, pprod)
            # curr_p was not used to build pprod, so it is an independent check
            if success && verifies_mod_p(new_qq_coeff, cnew_fq, F)
                pr.coeff_cands[i] = new_qq_coeff
                pr.is_stable[i] = true
                reg.n_unstable -= 1
                i += 1
                continue
            end
        end
        pr.mod_coeffs[i] = crt(ccurr, pprod, lift(ZZ, cnew_fq), curr_p)
        i += 1
    end
    reg.curr_ind += 1
    return ri
end

function update_registry!(reg::ModularRegistry{T}, p::T) where T
    push!(reg.pols, p)
    return length(reg.pols)
end
