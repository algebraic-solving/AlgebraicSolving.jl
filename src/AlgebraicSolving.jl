module AlgebraicSolving

greet() = print("AlgebraicSolving -- a package for algebraic solving based on msolve")

using msolve_jll

using AbstractAlgebra

import AbstractAlgebra:
    @attr,
    @attributes,
    @show_name,
    @show_special,
    addeq!,
    base_ring,
    canonical_unit,
    codomain,
    characteristic,
    data,
    degree,
    dim,
    domain,
    elem_type,
    evaluate,
    expressify,
    Field,
    FieldElem,
    force_coerce,
    force_op,
    gen,
    Generic,
    Generic.finish,
    Generic.MPolyBuildCtx,
    Generic.MPolyCoeffs,
    Generic.MPolyExponentVectors,
    Generic.push_term!,
    gens,
    get_attribute,
    get_attribute!,
    GF,
    Ideal,
    Map,
    map,
    MatElem,
    matrix,
    MatSpace,
    MPolyElem,
    MPolyRing,
    NCRing,
    NCRingElem,
    ngens,
    nvars,
    ordering,
    parent_type,
    PolyElem,
    PolynomialRing,
    PolyRing,
    Ring,
    RingElem,
    RingElement,
    set_attribute!,
    SetMap,
    symbols,
    total_degree

export PolynomialRing, MPolyRing, QQ, ZZ, FiniteField, GF, characteristic

include("interfaces/abstract-algebra.jl")

end # module AlgebraicSolving
