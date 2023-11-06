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
Pages = ["normal-forms.md"]
```

# Gröbner bases

## Introduction

AlgebraicSolving allows to compute normal forms of a polynomial resp. a finite
array of polynomials w.r.t. some given ideal over a finite field of
characteristic smaller $2^{31}$ w.r.t. the degree reverse lexicographical
monomial order.

**Note:** It therefore might first compute a Gröbner bases for the ideal.
## Functionality

```@docs
    normal_form(
        f::T,
        I::Ideal{T};
        nr_thrds::Int=1,
        info_level::Int=0
        ) where T <: MPolyRingElem
```

```@docs
    normal_form(
        F::Vector{T},
        I::Ideal{T};
        nr_thrds::Int=1,
        info_level::Int=0
        ) where T <: MPolyRingElem
```
