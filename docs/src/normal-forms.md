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

AlgebraicSolving allows to compute normal forms of a finite array of polynomials w.r.t.
some given ideal over finite
fields of characteristic smaller $2^{31}$ w.r.t. the degree reverse
lexicographical monomial order.

**Note:** It therefore might first compute a Gröbner bases for the ideal.
## Functionality

```@docs
    normal_form(
        F::Vector{T} where T <: MPolyRingElem,
        I::Ideal{T} where T <: MPolyRingElem;
        nr_thrds::Int=1,
        info_level::Int=0
        )
```
