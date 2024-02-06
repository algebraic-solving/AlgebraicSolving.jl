exit()
using Pkg; Pkg.activate("../../../")
using AlgebraicSolving
P, vrs = polynomial_ring(GF(65521), ["x$i" for i in 1:8]);
F = cyclic(P).gens;
