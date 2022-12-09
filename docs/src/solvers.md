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
Pages = ["solvers.md"]
```

# Algebraic Systems Solving

## Introduction

AlgebraicSolving allows to solve systems for input generators over finite
fields of characteristic smaller $2^{31}$ and over the rationals.

The underlying engine is provided by msolve.

## Functionality

```@docs
    rational_parametrization(
        I::Ideal{T} where T <: MPolyElem;
        initial_hts::Int=17,
        nr_thrds::Int=1,
        max_nr_pairs::Int=0,
        la_option::Int=2,
        info_level::Int=0,
        precision::Int=32
        )

    real_solutions(
        I::Ideal{T} where T <: MPolyElem;
        initial_hts::Int=17,
        nr_thrds::Int=1,
        max_nr_pairs::Int=0,
        la_option::Int=2,
        info_level::Int=0,
        precision::Int=32
        )
    rational_solutions(
        I::Ideal{T} where T <: MPolyElem;
        initial_hts::Int=17,
        nr_thrds::Int=1,
        max_nr_pairs::Int=0,
        la_option::Int=2,
        info_level::Int=0,
        precision::Int=32
        )
```

