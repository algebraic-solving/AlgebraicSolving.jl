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
Pages = ["decomposition.md"]
```

# Equidimensional Decomposition

## Introduction

AlgebraicSolving.jl allows to compute equidimensional decompositions
of polynomial ideals with defined over prime fields of characteristic
smaller than $2^{31}$ and over the rationals. 

More precisely, given a polynomial ideal $I$ it computes locally
closed sets $X_1, \dots, X_k$ s.t. $V(I)=\bigcup_{i=1}^{k} X_j$ and
such that each $X_j$ is equidimensional in the sense that the Zariski
closure of each $X_j$ is equidimensional.

The implemented algorithm is the one given in [this paper](https://arxiv.org/abs/2409.17785).

## Functionality

Each locally closed set in the output of the implemented decomposition
algorithm is of the form $V(F) \setminus V(g_1 \cdot \dots \cdot g_r)$
for a finite set of polynomials $F$ and polynomials $g_1,\dots, g_r$.

```@docs
    equations(
	    X::LocallyClosedSet
		)
```

```@docs
    inequations(
	    X::LocallyClosedSet
		)
```

```@docs
    Ideal(
		X::LocallyClosedSet
		)
```

```@docs
    equidimensional_decomposition(
	    I::Ideal{T};
	    info_level::Int=0
	    ) where {T <: MPolyRingElem}
```
