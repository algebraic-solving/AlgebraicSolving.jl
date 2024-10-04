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
Pages = ["dimension.md"]
```

# Krull dimension of an ideal

## Introduction

AlgebraicSolving allows to compute the Krull dimension for the ideal spanned
by given input generators over finite fields of characteristic smaller
$2^{31}$ and over the rationals.

The underlying engine is provided by msolve.

## Functionality

```@docs
    dimension(I::Ideal{T}) where T <: MPolyRingElem
```

