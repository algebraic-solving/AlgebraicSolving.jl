using Pkg; Pkg.activate("../../../")
using AlgebraicSolving
P, (x105, x104, x103, x102, x101, x30, x23, x22, x21, x13, x12, x11, x10) = polynomial_ring(GF(1301219), ["x105", "x104", "x103", "x102", "x101", "x30", "x23", "x22", "x21", "x13", "x12", "x11", "x10"],ordering=:degrevlex);
F = [x21-x12-x13, x22-x11-x13, x23-x11-x12, -x11^3-x12^3-x13^3+x30, x21*x22*x23-x10*x30, -3*x104*x11^2-x102-x103, -3*x104*x12^2-x101-x103, -3*x104*x13^2-x101-x103, x105*x22*x23+x101, x105*x22*x23+x101, x105*x21*x23+x102, x105*x21*x22+x103, x105*x21*x22+x103, -x10*x105+x104, -x105*x30+1];
sort!(F, by = p -> total_degree(p));
