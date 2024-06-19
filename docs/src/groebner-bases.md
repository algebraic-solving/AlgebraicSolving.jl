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
Pages = ["groebner-bases.md"]
```

# Gröbner bases

## Introduction

AlgebraicSolving allows to compute Gröbner bases for input generators over finite
fields of characteristic smaller $2^{31}$ w.r.t. the degree reverse
lexicographical monomial order.

At the moment different variants of Faugère's F4 Algorithm are implemented as
well as a signature based algorithm to compute Gröbner bases.

## Functionality

```@docs
    groebner_basis(
        I::Ideal{T} where T <: MPolyRingElem;
        initial_hts::Int=17,
        nr_thrds::Int=1,
        max_nr_pairs::Int=0,
        la_option::Int=2,
        eliminate::Int=0,
        intersect::Bool=true,
        complete_reduction::Bool=true,
        info_level::Int=0
        )
```

The engine supports the elimination of one block of variables considering the
product monomial ordering of two blocks, both ordered w.r.t. the degree
reverse lexicographical order. One can either directly add the number of
variables of the first block via the `eliminate` parameter in the
`groebner_basis` call. By using `intersect=false` it is possible to only 
use block ordering without intersecting. We have also implemented an alias 
for this call:

```@docs
    eliminate(
        I::Ideal{T} where T <: MPolyRingElem,
        eliminate::Int;
        intersect::Bool=true,
        initial_hts::Int=17,
        nr_thrds::Int=1,
        max_nr_pairs::Int=0,
        la_option::Int=2,
        complete_reduction::Bool=true,
        info_level::Int=0
        )
```

To compute signature Gröbner bases use

```@docs
    sig_groebner_basis(sys::Vector{T}; info_level::Int = 0, degbound::Int = 0) where {T <: MPolyRingElem}
```
