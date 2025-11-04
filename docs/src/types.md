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

Ideals can be constructed by giving an array of generators:

```@repl
using AlgebraicSolving
R, (x,y,z) = polynomial_ring(QQ, ["x", "y", "z"], internal_ordering=:degrevlex)
I = Ideal([x+y+1, y*z^2-13*y^2])
```

Ideals cache various data structures connected to ideals in order to make
computational algebra more effective. Assume that `T <: MPolyRingElem`, then
  * `gens::Vector{T}`: generators of the Ideal provded by the user;
  * `deg::Int`: degree of the ideal (`Nothing` if unknown);
  * `dim::Int`: Krull dimension of the Ideal (`Nothing` if unknown);
  * `gb::Dict{Int, Vector{T}}`: Gröbner bases of the Ideal where the key correspond to the number of variables eliminated (starting with the first);
  * `inter_sols::Vector{Vector{Vector{QQFieldElem}}}`: intervals with rational endpoints, containing the real solutions of the system (if `dim=0`);
  * `real_sols::Vector{Vector{QQFieldElem}}`: midpoints of the above intervals;
  * `rat_sols::Vector{Vector{QQFieldElem}}`: rational solutions of the system (if `dim=0`).
  * `rat_param::Union{RationalParametrization, RationalCurveParametrization}`: rational parametrization encoding the complex solution set of the Ideal (see below).
## Rational Parametrizations

Rational Parametrizations are special data structure for representing equidimensionnal
algebraic sets of dimension 0 and 1. They are typical outputs of algorithm for solving
polynomial systems encoding such algebraic sets thanks to their nice behaviour with
subsequent computations.

### Zero-dimensionnal parametrizations: finite set of points (see also [there](https://msolve.lip6.fr/downloads/msolve-tutorial.pdf#section.6))

A *zero-dimensional parametrization* $\mathscr{P}$ with
coefficients in a field $\mathbb{Q}$ consists of:
* polynomials $(w, w', v_1, \ldots, v_n)$ in $\mathbb{K}[t]$ where $t$ is
  a new variable and such that
    * $w$ is a monic square-free polynomial;
    * $w'=1$ when $\mathbb{K}$ is a prime field, $w'=\frac{dw}{dt}$ else;
    * $\deg(v_i) < \deg(w)$.
*  a linear form $l$  in the variables $x_1,\dotsc, x_n$,
such that
$$
l(\rho_1,\dotsc,\rho_n) = t \mod w.
$$

Such a data-structure encodes the following finite set of points
$$
\left\{\left (\frac{\rho_1(\vartheta)}{w'(\vartheta)},
\ldots, \frac{\rho_n(\vartheta)}{w'(\vartheta)} \right) \middle|\, w(\vartheta) = 0\right\}.
$$
According to this definition, the roots of $w$ are exactly the
values taken by $l$ on this set.

The type `RationalParametrization` therefore caches the following attributes:
  * `vars::Vector{Symbol}`: variables used for the parametrization the last one playing the role of the above $t$
  (hence, maybe with one more variable than the ones given as input);
  * `cfs_lf::Vector{ZZRingElem}`: coefficients on the linear form $l$ (when variables are added);
  * `elim::QQPolyRingElem`: elimination polynomial $w$ ;
  * `denom::QQPolyRingElem`: denominator polynomial (usually $w'$);
  * `param::Vector{QQPolyRingElem}`: numerators $\rho_i$'s.

See the documentation of the [rational_parametrization](@ref) function for further details.

### One-dimensionnal parametrizations: algebraic curves

A *one-dimensional parametrization* $\mathscr{C}$ with
coefficients in a field $\mathbb{C}$ consists of:
* polynomials $(w, w', v_1, \ldots, v_n)$ in $\mathbb{K}[t,s]$ where $t,s$ are
  new variables and such that
    * $w$ is a monic square-free polynomial and $w'=\frac{\partial w}{\partial t}$;
    * $\deg(v_i) < \deg(w)$.
*  two linear forms $(l,l')$  in the variables $x_1,\dotsc, x_n$,
such that
$$
l(\rho_1,\dotsc,\rho_n) = t\, w'\mod w \quad \text{and} \quad
l'(\rho_1,\dotsc,\rho_n) = s\,w' \mod w
$$

Such a data-structure encodes the curve defined as the Zariski closure of the following set
$$
\left\{\left (\frac{\rho_1(\vartheta,\eta)}{w'(\vartheta,\eta)},
\ldots, \frac{\rho_n(\vartheta,\eta)}{w'(\vartheta,\eta)} \right)
\middle|\, w(\vartheta,\eta) = 0,\, w'(\vartheta, \eta) \neq 0\right\}.
$$
According to this definition, the roots of $w$ are exactly the
values taken by $l$ on this set.

The type `RationalCurveParametrization` therefore caches the following attributes:
  * `vars::Vector{Symbol}`: variables used for the parametrization the last two ones playing the role of the above $(t,s)$
  (hence, maybe with up to two more variable than the ones given as input);;
  * `cfs_lfs::Vector{Vector{ZZRingElem}}`: coefficients on the linear form $l$ (when variables are added);
  * `elim::QQMPolyRingElem`: elimination polynomial $w$ ;
  * `denom::QQMPolyRingElem`: denominator polynomial (usually $w'$);
  * `param::Vector{QQMPolyRingElem}`: numerators $\rho_i$'s.

See the documentation of the [rational_curve_parametrization](@ref) function for further details.

## Roadmaps

Roadmaps are algebraic curves capturing the connectivity properties of an algebraic
set containing it. They are computed as union of several curves, organized in a tree
due to the recursive nature of the roadmap algorithms.

```julia
mutable struct Roadmap
    initial_ideal::Ideal{QQMPolyRingElem}
    root::RMnode
end
```

Hence, it is composed of a main node, containing the equations of the initial algebraic set and a reference to the root node of the tree.

```julia
mutable struct RMnode
    polar_eqs::Vector{QQMPolyRingElem}
    base_pt::Vector{QQFieldElem}
    children::Vector{RMnode}
end
```

Each roadmap node (including the root) correspond to a curve component of the roadmap. These are defined in fibers of the initial variety, and are arranged w.r.t. to the specialized variables they share. In particular, all child nodes (referenced in `children`) of a node will have the same specialized variable, plus a new distinct one.

Each of these components is defined as the polar variety (whose equations are in `polar_eqs`) over a fiber of the initial algebraic set, which are defined by specializing the first variables to the `base_pt` of its parents, and of itsef (in this order).


See the documentation of the [roadmap](@ref) function for further details.


```
TODO : un dessin ? les fonctions all_eqs, all_base_pts et nb_nodes (utiliser la doc du code?)
```