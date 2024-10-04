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

AlgebraicSolving allows to compute equidimensional decompositions of polynomial
ideals. This is to be understood in a geometric sense, i.e. given a polynomial
ideal $I$ it computes ideals $I_1,\dots,I_k$ s.t. $V(I)=\bigcup_{i=1}^{k} V(I_j)$
and such that each $V(I_j)$ is equidimensional.

The implemented algorithm is the one given in [this paper](https://arxiv.org/abs/2409.17785).

## Functionality

```@docs
    equidimensional_decomposition(
	    I::Ideal{T}, 
		info_level::Int=0
		) where {T <: MPolyRingElem}
```
