# Getting Started

AlgebraicSolving is a computer algebra package for the Julia programming 
language, maintained by Christian Eder, Mohab Safey El Din, Rafael Mohr, Rémi Prébet.

- <https://github.com/algebraic-solving/AlgebraicSolving.jl> (Source code)

The features of AlgebraicSolving include algorithms for computing
Gröbner bases over finite fields and for computing real solutions.
The main workhorse of AlgebraicSolving is the [msolve
library](https://msolve.lip6.fr/) .

## Installation

AlgebraicSolving requires Julia 1.10 or higher. Please see
<https://julialang.org/downloads/> for instructions on
how to obtain julia for your system.

To install AlgebraicSolving, simply type

```
julia> using Pkg; Pkg.add("AlgebraicSolving")
```
in the Julia prompt.
