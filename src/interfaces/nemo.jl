
@doc Markdown.doc"""
    _convert_to_msolve(
            F::Vector{T}) where T <: MPolyRingElem

Convert a vector of polynomials to input data for msolve.

**Note**: This is an internal function.
"""
function _convert_to_msolve(
        F::Vector{T}) where T <: MPolyRingElem

    R = first(F).parent
    
    nr_vars    = nvars(R)
    nr_gens    = length(F)
    lens       = Int32[F[i].length for i in 1:nr_gens]
    nr_terms   = sum(lens)
    field_char = characteristic(R)

    if field_char > 2^31 || degree(base_ring(R)) != 1
        error("At the moment we only support prime fields up to prime characteristic < 2^31.")
    end
    # get coefficients
    if field_char == 0
        cfs = BigInt[]
    else
        cfs = Int32[]
    end
    if field_char == 0
        for i in 1:nr_gens
            for cf in coefficients(F[i])
                push!(cfs, BigInt(numerator(cf)))
                push!(cfs, BigInt(denominator(cf)))
            end
        end
    else
        for i in 1:nr_gens
            for cf in coefficients(F[i])
                push!(cfs, Int32(data(prime_field(base_ring(R))(cf))))
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

    return lens, cfs, exps
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
function _convert_finite_field_gb_to_abstract_algebra(
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

    basis = []
    #= basis = Vector{MPolyRingElem}(undef, bld) =#

    len   = 0

    if eliminate > 0
        z = zeros(Int, eliminate)
    end
    for i in 1:nr_gens
        #= check if element is part of the eliminated basis =#
        if eliminate > 0
            cmp = convert(Vector{Int}, bexp[(len)*nr_vars+1:(len+1)*nr_vars])
            if cmp[1:eliminate] > z
                continue
            end
        end
        g  = MPolyBuildCtx(R)
        for j in 1:blen[i]
            push_term!(g, CR(bcf[len+j]),
                       convert(Vector{Int}, bexp[(len+j-1)*nr_vars+1:(len+j)*nr_vars]))
        end
        len +=  blen[i]
        push!(basis, g.poly)
    end

    return basis
end
