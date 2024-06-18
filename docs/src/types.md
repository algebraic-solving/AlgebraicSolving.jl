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
Pages = ["types.md"]
```

# Data Types

## Introduction

AlgebraicSolving handles ideals in multivariate polynomial rings over a prime 
field of characteristic $0$ or $p$ where $p$ is a prime number $<2^{31}$.

## Polynomial Rings

We use [Nemo](https://www.nemocas.org/)'s multivariate polynomial 
ring structures:

```@repl
using AlgebraicSolving
R, (x,y,z) = polynomial_ring(QQ, ["x", "y", "z"], internal_ordering=:degrevlex)
```
The above example defines a multivariate polynomial ring in three variables `x`, 
`y`, and `z` over the rationals using the dgree reverse lexicographical ordering 
for printing polynomials in the following. One can also define polynomial rings 
over finite fields:

```@repl
using AlgebraicSolving
R, (x,y,z) = polynomial_ring(GF(101), ["x", "y", "z"], internal_ordering=:degrevlex)
```

## Ideals

Ideals can be constructed by giving an array of generators. Ideals cache varies 
data structures connected to ideals in order to make computational algebra more 
effective:

```@repl
using AlgebraicSolving
R, (x,y,z) = polynomial_ring(QQ, ["x", "y", "z"], internal_ordering=:degrevlex)
I = Ideal([x+y+1, y*z^2-13*y^2])
```

