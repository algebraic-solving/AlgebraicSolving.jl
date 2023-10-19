@doc Markdown.doc"""
    cyclic(R::MPolyRing)

Return the Cyclic ideal in the variables of `R`.

# Example
```jldoctest
julia> using AlgebraicSolving

julia> R, vars = polynomial_ring(QQ, ["x$i" for i in 1:4])
(Multivariate polynomial ring in 4 variables over QQ, QQMPolyRingElem[x1, x2, x3, x4])

julia> cyclic(R)
QQMPolyRingElem[x1 + x2 + x3 + x4, x1*x2 + x1*x4 + x2*x3 + x3*x4, x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4, x1*x2*x3*x4 - 1]
```
"""
function cyclic(R::MPolyRing)
    vars = gens(R)
    n = length(vars)
    pols = [sum(prod(vars[j%n+1] for j in k:k+i) for k in 1:n) for i in 0:n-2]
    push!(pols, prod(vars[i] for i in 1:n)-1)
    return Ideal(pols)
end
