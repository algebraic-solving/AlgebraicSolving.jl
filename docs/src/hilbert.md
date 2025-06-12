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
Pages = ["hilbert.md"]
```

# Hilbert series of an ideal

## Introduction

AlgebraicSolving allows to compute the Hilbert series for the ideal spanned
by given input generators over finite fields of characteristic smaller
$2^{31}$ and over the rationals.

The underlying engine is provided by msolve.

## Functionality

```@docs
    hilbert_series(I::Ideal{T}) where T <: MPolyRingElem
```

In addition, from the same input, AlgebraicSolving can also compute the dimension and degree of the ideal, as well as the Hilbert polynomial and index of regularity.

```@docs
    hilbert_dimension(I::Ideal{T}) where T <: MPolyRingElem

    hilbert_degree(I::Ideal{T}) where T <: MPolyRingElem

    hilbert_polynomial(I::Ideal{T}) where T <: MPolyRingElem
```