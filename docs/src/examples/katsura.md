```@meta
CurrentModule = AlgebraicSolving
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

