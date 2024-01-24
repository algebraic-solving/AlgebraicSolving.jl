using Pkg; Pkg.activate("../../../")
using AlgebraicSolving
P, (x, y, z) = polynomial_ring(GF(65521), ["x", "y", "z"]);
F = [x*y, x*z, y*z];
sort!(F, by = p -> total_degree(p));
