```@meta
CurrentModule = AlgebraicSolving
DocTestSetup = quote
  using AlgebraicSolving
end
```

```@setup algebraicsolving
using AlgebraicSolving
```

```@contents
Pages = ["connectivity.md"]
```

# Algorithms for connectivity queries on smooth bounded algebraic sets

## Introduction

AlgebraicSolving allows to compute a roadmap for the real trace of the zero-set
of the ideal spanned by given input generators over the rationals.

It assumes that the underlying algebraic set is **smooth**, and its real trace is **bounded**.

The underlying engine is provided by msolve.

## Functionality

```@docs
    roadmap(
      I::Ideal{P} where P <: QQMPolyRingElem;
      C::Vector{Vector{Vector{QQFieldElem}}}=Vector{Vector{QQFieldElem}}[],
      info_level::Int=0,
      checks::Bool=false
    )
```

In addition, AlgebraicSolving can compute equations definition critical loci of polynomial maps over the given algebraic set.

```@docs
   computepolar(
        J::Union{Vector{Int},UnitRange{Int}},
        V::Ideal{P};
        phi::Vector{P} = P[],
        dimproj = length(J)-1,
        only_mins = false
    ) where (P <: MPolyRingElem)
```