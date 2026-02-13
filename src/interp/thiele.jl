struct ThieleError <: Exception
    msg::String
end

Base.show(io::IO, e::ThieleError) = print(io, "ThieleError: ", e.msg)

struct ThieleDivideError <: Exception end

function thiele(
    R::QQPolyRing,
    bb::Function;
    retry::Int=10,
    show_progress::Bool=false,
    offset::Int=1
)::FracFieldElem{QQPolyRingElem}
    x = QQFieldElem[]
    y = QQFieldElem[]
    # We denote by ρ(i, n) := ρₙ(xᵢ, xᵢ₊₁, ..., xᵢ₊ₙ)
    # the reciprocal difference of order n of the points (xᵢ, xᵢ₊₁, ..., xᵢ₊ₙ).
    # Following https://rosettacode.org/wiki/Thiele%27s_interpolation_formula#C,
    # it can be computed recursively as
    ρ_cache = Dict{Tuple{Int,Int},QQFieldElem}()
    function ρ(i::Int, n::Int)
        if n == -2 || n == -1
            return QQ(0)
        elseif n == 0
            return y[i]
        elseif haskey(ρ_cache, (i, n))
            # We cache the computed values to avoid redundant computations
            return ρ_cache[(i, n)]
        end
        den = ρ(i, n - 1) - ρ(i + 1, n - 1)
        if den == 0
            # ρ(i, n) is not well-defined.
            # This happens when the new point depends on previous points.
            # It is very likely that the interpolating formula
            # does not require a new reciprocal difference.
            # We stop here and return the formula with n-1 points.
            throw(ThieleDivideError())
        end
        val = (x[i] - x[i+n]) / den + ρ(i + 1, n - 2)
        ρ_cache[(i, n)] = val
        val
    end
    # We keep adding points until we have enough points to compute
    # the interpolating formula.
    prog = ProgressBar(desc="Univariate rational interpolation", offset=offset, enabled=show_progress)
    update!(prog, 0)
    n = 1
    xᵢ = QQ(n)
    # We allow the black box function to fail for some points.
    # But, if it fails for too many consecutive points, we assume that
    # the function is not well-defined and stop the interpolation.
    failures = 0
    while true
        xᵢ = xᵢ + 1
        yᵢ = try
            bb(xᵢ)
        catch e
            if isa(e, DivideError)
                failures += 1
                if failures >= retry
                    throw(ThieleError("Too many consecutive failures evaluating the black box function."))
                end
                continue
            else
                rethrow(e)
            end
        end
        failures = 0
        push!(y, yᵢ)
        push!(x, xᵢ)
        try
            ρ(1, n - 1)
        catch e
            if isa(e, ThieleDivideError)
                break
            else
                rethrow(e)
            end
        end
        next!(prog)
        n += 1
    end
    finish!(prog)
    # Now that we have collected enough points, we can compute
    # the interpolating formula.
    function F(n::Int, N::Int)
        if n == N
            # Note that in the link, F_N(x) is defined incorrectly.
            return fraction_field(R)(ρ(1, N) - ρ(1, N - 2))
        else
            return ρ(1, n) - ρ(1, n - 2) + (gen(R) - x[n+1]) / F(n + 1, N)
        end
    end
    f = try
        F(0, n - 2)
    catch e
        if isa(e, ThieleDivideError) || isa(e, DivideError)
            throw(ThieleError("Division by zero encountered while computing the interpolating formula."))
        else
            rethrow(e)
        end
    end
    # We verify that the interpolating formula is correct
    # for an additional point.
    xᵢ = xᵢ + 1
    if evaluate(f, xᵢ) != bb(xᵢ)
        throw(ThieleError("Verification failed."))
    end
    return f
end
