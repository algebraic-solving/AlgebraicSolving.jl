```@meta
CurrentModule = AlgebraicSolving
```

# Data Types

## Introduction

AlgebraicSolving handles ideals in multivariate polynomial rings over a prime 
field of characteristic $0$ or $p$ where $p$ is a prime number $<2^{31}$.

## Polynomial Rings

We use [Nemo](https://www.nemocas.org/index.html)'s multivariate polynomial 
ring structures:

```@repl
using AlgebraicSolving
R, (x,y,z) = PolynomialRing(QQ, ["x", "y", "z"], ordering=:degrevlex)
```
The above example defines a multivariate polynomial ring in three variables `x`, 
`y`, and `z` over the rationals using the dgree reverse lexicographical ordering 
for printing polynomials in the following. One can also define polynomial rings 
over finite fields:

```@repl
using AlgebraicSolving
R, (x,y,z) = PolynomialRing(GF(101), ["x", "y", "z"], ordering=:degrevlex)
```

## Ideals

Ideals can be constructed by giving an array of generators. Ideals cache varies 
data structures connected to ideals in order to make computational algebra more 
effective:

```@repl
using AlgebraicSolving
R, (x,y,z) = PolynomialRing(QQ, ["x", "y", "z"], ordering=:degrevlex)
I = Ideal([x+y+1, y*z^2-13*y^2])
```

