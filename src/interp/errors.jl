@doc Markdown.doc"""
    UnluckyPointError(msg::String)

Exception thrown by a black box function to signal that it cannot be evaluated
meaningfully at the given point, e.g. because the point lies on a proper
algebraic subset of the parameter space where the specialized object has a
different structure than at a generic point.

The interpolation routines (`thiele`, `cuyt_lee`) treat such a point as unlucky:
they skip it and evaluate the black box function elsewhere, within their retry
limits. Any other exception thrown by a black box function is propagated.
"""
struct UnluckyPointError <: Exception
    msg::String
end

Base.show(io::IO, e::UnluckyPointError) = print(io, "UnluckyPointError: ", e.msg)
