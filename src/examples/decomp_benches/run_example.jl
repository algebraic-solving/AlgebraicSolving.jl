Fhom = homogenize(F);
S = parent(first(Fhom));
x, y, z = AlgebraicSolving.gens(S)[1:3]
tst = [x*y, x*z, y*z];
sig_decomp(tst, info_level = 1);
println("----------")
tm = @elapsed bla = sig_decomp(Fhom, info_level=1);
println("time $(tm), $(length(bla)) components")
