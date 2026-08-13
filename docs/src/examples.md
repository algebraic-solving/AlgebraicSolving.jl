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
Pages = ["examples.md"]
```

# Examples

Here we include some well-known example multivariate polynomial systems.

## Katsura-n

These systems appeared in a problem of magnetism in physics.
For a given $n$ `katsura(n)` has $2^n$ solutions and is defined in a
polynomial ring with $n+1$ variables. For a given polynomial ring `R`
with $n$ variables `katsura(R)` defines the corresponding system with
$2^{n-1}$ solutions.

### Functionality

```@docs
    katsura(n::Int)
    katsura(R::MPolyRing)
```

## Eco-n

These systems appeared in an economics modeling problem described by Alexander Morgan.
For a given $n \geq 1$ `eco(n)` is defined in a polynomial ring with $n$ variables.
For $n \geq 2$, it has $2^{n-2}$ solutions, while for $n = 1$ it degenerates to the
unit ideal. For a given polynomial ring `R` with $n$ variables, `eco(R)` defines the
corresponding system.

### Functionality

```@docs
    eco(n::Int)
    eco(R::MPolyRing)
```
