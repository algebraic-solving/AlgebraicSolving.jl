```@meta
CurrentModule = AlgebraicSolving
```

# Gröbner bases

## Introduction

AlgebraicSolving allows to compute Gröbner bases for input generators over finite
fields of characteristic smaller $2^{31}$ w.r.t. the degree reverse
lexicographical monomial order.

At the moment different variants of Faugère's F4 Algorithm are implemented.

## Functionality

```@docs
    groebner_basis(
        I::Ideal{T} where T <: MPolyElem;
        initial_hts::Int=17,
        nr_thrds::Int=1,
        max_nr_pairs::Int=0,
        la_option::Int=2,
        eliminate::Int=0,
        complete_reduction::Bool=true,
        info_level::Int=0
        )
```

The engine supports the elimination of one block of variables considering the
product monomial ordering of two blocks, both ordered w.r.t. the degree
reverse lexicographical order. One can either directly add the number of
variables of the first block via the `eliminate` parameter in the
`groebner_basis` call. We have also implemented an alias for this call:

%% ```@docs
%% function eliminate(
%%         I::Ideal{T} where T <: MPolyElem,
%%         eliminate::Int,
%%         initial_hts::Int=17,
%%         nr_thrds::Int=1,
%%         max_nr_pairs::Int=0,
%%         la_option::Int=2,
%%         complete_reduction::Bool=true,
%%         info_level::Int=0
%%         )
%% ```
