function _resize_ff!(a::FqMPolyRingElem, n::Int)
    if a.data isa Generic.MPoly
        fit!(a, n)
    elseif a.data isa fqPolyRepMPolyRingElem
        ccall((:fq_nmod_mpoly_resize, Nemo.libflint), Nothing,
            (Ref{fqPolyRepMPolyRingElem}, Int, Ref{fqPolyRepMPolyRing}), a.data,
            n, a.parent.data)
    else
        @assert a.data isa Nemo.fpMPolyRingElem
        ccall((:nmod_mpoly_resize, Nemo.libflint), Nothing,
            (Ref{fpMPolyRingElem}, Int, Ref{fpMPolyRing}), a.data, n,
            a.parent.data)
     end
end

function _resize_ff!(a::fpMPolyRingElem, n::Int)
    ccall((:nmod_mpoly_resize, Nemo.libflint), Cvoid, (Ref{fpMPolyRingElem}, Int, Ref{fpMPolyRing}), a, n, parent(a))
end

function _resize_qq!(a::QQMPolyRingElem, n::Int)
    ccall((:fmpq_mpoly_resize, Nemo.libflint), Cvoid, (Ref{QQMPolyRingElem}, Int, Ref{QQMPolyRing}), a, n, parent(a))
end

@doc Markdown.doc"""
    _convert_to_msolve(
            F::Vector{T}) where T <: MPolyRingElem

Convert a vector of polynomials to input data for msolve.

**Note**: This is an internal function.
"""
function _convert_to_msolve(
        F::Vector{T}) where T <: MPolyRingElem

    R = parent(first(F))
    
    nr_vars    = nvars(R)
    nr_gens    = length(F)
    lens       = Int32[length(F[i]) for i in 1:nr_gens]
    nr_terms   = sum(lens)
    field_char = characteristic(R)

    if field_char > 2^31 || degree(base_ring(R)) != 1
        error("At the moment we only support prime fields up to prime characteristic < 2^31.")
    end

    lens = Int32[]
    # get coefficients
    if field_char == 0
        cfs = BigInt[]
    else
        cfs = Int32[]
    end
    if field_char == 0
        for i in 1:nr_gens
            if F[i] != R(0)
                for cf in coefficients(F[i])
                    push!(cfs, BigInt(numerator(cf)))
                    push!(cfs, BigInt(denominator(cf)))
                end
                push!(lens, length(F[i]))
            end
        end
    else
        for i in 1:nr_gens
            if F[i] != R(0)
                for cf in coefficients(F[i])
                    push!(cfs, Int32(lift(Nemo.ZZ, cf)))
                end
                push!(lens, length(F[i]))
            end
        end
    end

    # get exponent vectors
    exps  = Int32[]
    for i in 1:nr_gens
        for ev in exponent_vectors(F[i])
            append!(exps, convert(Vector{Int32}, ev))
        end
    end

    return lens, cfs, exps, length(lens)
end

@doc Markdown.doc"""
    _convert_finite_field_gb_to_abstract_algebra(
        bld::Int32,
        blen::Vector{Int32},
        bcf::Vector{Int32},
        bexp::Vector{Int32},
        R::MPolyRing,
        eliminate::Int=0
        )

Converts a finite Gröbner basis computed internally via msolve
to a vector of polynomials.

**Note**: This is an internal function.
"""
function _convert_finite_field_array_to_abstract_algebra(
        bld::Int32,
        blen::Vector{Int32},
        bcf::Vector{Int32},
        bexp::Vector{Int32},
        R::MPolyRing,
        eliminate::Int=0
        )

    if characteristic(R) == 0
        error("At the moment we only support finite fields.")
    end

    nr_gens = bld
    nr_vars = nvars(R)
    CR      = coefficient_ring(R)

    len = 0
    ctr = 0

    if eliminate > 0
        z = zeros(Int, eliminate)
    end
    g = [zero(R) for j in 1:nr_gens]
    for i in 1:nr_gens
        #= check if element is part of the eliminated basis =#
        if eliminate > 0
            cmp = convert(Vector{Int}, bexp[(len)*nr_vars+1:(len+1)*nr_vars])
            if cmp[1:eliminate] > z
                continue
            end
        end
        if bcf[len+1] == 0
            g[i] = R(0)
        else
            _resize_ff!(g[i], Int(blen[i]))
            for j in 1:blen[i]
                Nemo.setcoeff!(g[i], j, CR(bcf[len+j]))
            end
            for j in 1:blen[i]
                Nemo.set_exponent_vector!(g[i], j, convert(Vector{Int}, bexp[(len+j-1)*nr_vars+1:(len+j)*nr_vars]))
            end
        end
        ctr += 1
        len +=  blen[i]
    end
    basis = g[1:ctr]
    sort_terms!.(basis)
    return basis
end

function _convert_rational_array_to_abstract_algebra(
        bld::Int32,
        blen::Vector{Int32},
        bcf::Vector{QQFieldElem},
        bexp::Vector{Int32},
        R::MPolyRing,
        normalize::Bool=false,
        eliminate::Int=0
        )

    #  Note: Over QQ msolve already returns the eliminated basis
    #  only in order to save memory due to possible huge
    #  coefficient sizes.
    if characteristic(R) != 0
        error("We assume QQ as base field here.")
    end

    nr_gens = bld
    nr_vars = nvars(R)
    CR      = coefficient_ring(R)

    len   = 0

    g = [zero(R) for j in 1:nr_gens]
    for i in 1:nr_gens
        if bcf[len+1] == 0
            g[i] = R(0)
        else
            _resize_qq!(g[i], Int(blen[i]))
            lc = bcf[len+1]

            if normalize && lc != 1
                for j in 1:blen[i]
                    Nemo.setcoeff!(g[i], j, bcf[len+j]/lc)
                end
            else
                for j in 1:blen[i]
                    Nemo.setcoeff!(g[i], j, bcf[len+j])
                end
            end
            for j in 1:blen[i]
                Nemo.set_exponent_vector!(g[i], j, convert(Vector{Int}, bexp[(len+j-1)*nr_vars+1:(len+j)*nr_vars]))
            end
        end
        len +=  blen[i]
    end

    sort_terms!.(g)
    return g
end
