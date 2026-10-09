# Benchmark for `equidimensional_decomposition` on the "soscritpoints_deg_2_nvars_5_s_2" system.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AlgebraicSolving
using FileWatching.Pidfile: mkpidlock

const NAME = "soscritpoints_deg_2_nvars_5_s_2"

R, (x1, x2, x3, x4, x5) = polynomial_ring(QQ, ["x1", "x2", "x3", "x4", "x5"])

F = [
    200*x1^2*x2 + 1640*x1^2*x3 + 5348*x1^2*x5 - 960*x1*x2^2 - 5248*x1*x2*x3 + 4544*x1*x2*x5 - 1660*x1*x3*x4 - 8178*x1*x3*x5 - 180*x1*x4*x5 - 9024*x1*x5 + 1024*x2^3 + 5312*x2*x3*x4 + 576*x2*x4*x5 + 4418*x2*x5^2 - 8460*x4*x5^2 + 4042*x5^3,
    -12528*x1^3 + 1640*x1^2*x2 + 28586*x1^2*x3 - 11644*x1^2*x5 + 16704*x1^2 - 2624*x1*x2^2 - 1660*x1*x2*x4 - 8178*x1*x2*x5 - 27224*x1*x3*x4 + 25970*x1*x4*x5 - 7482*x1*x5^2 + 2656*x2^2*x4 + 13778*x3*x4^2 + 1494*x4^2*x5,
    -12960*x1^2*x5 - 1660*x1*x2*x3 - 180*x1*x2*x5 - 13612*x1*x3^2 + 25970*x1*x3*x5 + 1278*x1*x5^2 + 17280*x1*x5 + 2656*x2^2*x3 + 288*x2^2*x5 - 8460*x2*x5^2 + 13778*x3^2*x4 + 2988*x3*x4*x5 + 16362*x4*x5^2 - 7740*x5^3,
    5348*x1^2*x2 - 11644*x1^2*x3 - 12960*x1^2*x4 + 22466*x1^2*x5 + 2272*x1*x2^2 - 8178*x1*x2*x3 - 180*x1*x2*x4 - 9024*x1*x2 + 25970*x1*x3*x4 - 14964*x1*x3*x5 + 2556*x1*x4*x5 + 17280*x1*x4 - 16512*x1*x5 + 288*x2^2*x4 + 4418*x2^2*x5 - 16920*x2*x4*x5 + 12126*x2*x5^2 + 1494*x3*x4^2 + 16362*x4^2*x5 - 23220*x4*x5^2 + 7396*x5^3,
    5184*x1^4 - 12528*x1^3*x3 - 13824*x1^3 + 100*x1^2*x2^2 + 1640*x1^2*x2*x3 + 5348*x1^2*x2*x5 + 14293*x1^2*x3^2 - 11644*x1^2*x3*x5 + 16704*x1^2*x3 - 12960*x1^2*x4*x5 + 11233*x1^2*x5^2 + 9216*x1^2 - 320*x1*x2^3 - 2624*x1*x2^2*x3 + 2272*x1*x2^2*x5 - 1660*x1*x2*x3*x4 - 8178*x1*x2*x3*x5 - 180*x1*x2*x4*x5 - 9024*x1*x2*x5 - 13612*x1*x3^2*x4 + 25970*x1*x3*x4*x5 - 7482*x1*x3*x5^2 + 1278*x1*x4*x5^2 + 17280*x1*x4*x5 - 8256*x1*x5^2 + 256*x2^4 + 2656*x2^2*x3*x4 + 288*x2^2*x4*x5 + 2209*x2^2*x5^2 - 8460*x2*x4*x5^2 + 4042*x2*x5^3 + 6889*x3^2*x4^2 + 1494*x3*x4^2*x5 + 8181*x4^2*x5^2 - 7740*x4*x5^3 + 1849*x5^4
]

# Reducing the rational system is far cheaper than building the same
# polynomials again over the prime field.
prime = Int32(AlgebraicSolving.Nemo.rand_bits_prime(ZZ, 31))
Rp, _ = polynomial_ring(GF(prime), ["x1", "x2", "x3", "x4", "x5"])
Fp = [AlgebraicSolving.reduce_mod_p(f, Rp) for f in F]

# Decompose a tiny system first so that the timings below measure the
# computation rather than compilation. The solver specialises on the number of
# variables, so the warm up uses as many of them as the benchmark itself.
function warmup(K)
    _, w = polynomial_ring(K, ["w$(i)" for i in 1:5])
    equidimensional_decomposition(Ideal([w[1]*w[2], w[1]*w[3], w[2]*w[3]]))
    return nothing
end

warmup(QQ)
warmup(GF(prime))

# A short first run can still be dominated by one off compilation that the warm
# up does not cover, so measure a second one and keep that timing. Past this
# threshold the computation dwarfs the overhead and the first timing stands.
println("### $(NAME): equidimensional_decomposition over QQ ###")
time_qq = @elapsed equidimensional_decomposition(Ideal(F), info_level = 1)
if time_qq < 120
    println("### $(NAME): second run over QQ ###")
    time_qq = @elapsed equidimensional_decomposition(Ideal(F), info_level = 1)
end

println("### $(NAME): equidimensional_decomposition over GF($(prime)) ###")
time_gf = @elapsed equidimensional_decomposition(Ideal(Fp), info_level = 1)
if time_gf < 120
    println("### $(NAME): second run over GF($(prime)) ###")
    time_gf = @elapsed equidimensional_decomposition(Ideal(Fp), info_level = 1)
end

ratio = time_qq / time_gf

println()
println("### $(NAME) timings ###")
println("QQ:          $(time_qq) s")
println("GF($(prime)): $(time_gf) s")
println("QQ / GF:     $(ratio)")

# Collect the timings of every benchmark in this folder in one table. Rows are
# appended, so re-running a benchmark adds a row rather than replacing one.
# Concurrent runs would otherwise lose and interleave rows, hence the lock.
csv = joinpath(@__DIR__, "timings.csv")
mkpidlock(csv * ".lock") do
    write_header = !isfile(csv)
    open(csv, "a") do io
        write_header && println(io, "system,nvars,neqns,time_qq,time_gf,ratio_qq_gf")
        println(io, join((NAME, 5, 5,
                          round(time_qq, digits = 6),
                          round(time_gf, digits = 6),
                          round(ratio, digits = 4)), ","))
    end
end
println("appended to $(csv)")
