```@meta
CurrentModule = AlgebraicSolving
```

# Gröbner bases

## Introduction

AbstractAlgebra allows to compute Gröbner bases for input generators over finite
fields of characteristic smaller $2^{31}$.

At the moment different variants of Faugère's F4 Algorithm are implemented.

## Functionality

```@docs
    groebner_basis(
        F::Vector{T} where T <: MPolyElem;
        initial_hts::Int=17,
        nr_thrds::Int=1,
        max_nr_pairs::Int=0,
        la_option::Int=2,
        eliminate::Int=0,
        complete_reduction::Bool=true,
        info_level::Int=0
        )
```

