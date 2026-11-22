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

## Chandra-n

These systems come from a discretization of Chandrasekhar's H-equation, which
appears in radiative transfer. For a given $n \geq 1$ `chandra(n)` is defined in a
polynomial ring with $n$ variables and has $2^{n-1}$ solutions. The parameter `c`
is usually chosen as $0.51234$. For a given polynomial ring `R` with $n$ variables,
`chandra(R)` defines the corresponding system.

### Functionality

```@docs
    chandra(n::Int, characteristic::Int=0)
    chandra(R::MPolyRing)
```
