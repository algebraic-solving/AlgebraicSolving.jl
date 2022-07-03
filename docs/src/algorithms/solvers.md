```@meta
CurrentModule = AlgebraicSolving
```

# Algebraic Systems Solving

## Introduction

AbstractAlgebra allows to solve systems for input generators over finite
fields of characteristic smaller $2^{31}$ and over the rationals.

The underlying engine is provided by msolve.

## Functionality

```@docs
    rational_parametrization(
        F::Vector{T} where T <: MPolyElem;
        initial_hts::Int=17,
        nr_thrds::Int=1,
        max_nr_pairs::Int=0,
        la_option::Int=2,
        info_level::Int=0,
        precision::Int=32
        )

    real_solutions(
        F::Vector{T} where T <: MPolyElem;
        initial_hts::Int=17,
        nr_thrds::Int=1,
        max_nr_pairs::Int=0,
        la_option::Int=2,
        info_level::Int=0,
        precision::Int=32
        )
```

