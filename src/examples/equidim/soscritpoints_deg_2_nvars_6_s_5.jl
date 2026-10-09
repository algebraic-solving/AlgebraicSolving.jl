# Benchmark for `equidimensional_decomposition` on the "soscritpoints_deg_2_nvars_6_s_5" system.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AlgebraicSolving
using FileWatching.Pidfile: mkpidlock

const NAME = "soscritpoints_deg_2_nvars_6_s_5"

R, (x1, x2, x3, x4, x5, x6) = polynomial_ring(QQ, ["x1", "x2", "x3", "x4", "x5", "x6"])

F = [
    -1156*x1*x2*x5 - 1054*x1*x3*x5 + 15080*x2^3 + 3162*x2^2*x3 + 19116*x2^2*x4 + 7972*x2*x3^2 - 3740*x2*x3*x5 - 16280*x2*x3 + 2296*x2*x4^2 + 12272*x2*x4 + 22002*x2*x5^2 - 3700*x2*x5 + 5192*x2*x6^2 - 2108*x2*x6 + 3636*x2 + 9130*x3^2*x5 - 3224*x3*x4^2 - 5998*x3*x4*x5 - 3080*x3*x4*x6 - 94*x3*x5^2 - 28412*x3*x5 - 1922*x3*x6 + 14592*x4^2*x5 + 5616*x4^2 + 20450*x4*x5^2 + 952*x4*x5*x6 + 5624*x4*x5 + 2376*x4*x6^2 + 4144*x4*x6 - 3348*x4 - 15936*x5^2,
    72*x1^2*x3 - 636*x1^2 - 1054*x1*x2*x5 - 1248*x1*x3*x6 + 480*x1*x4^2 - 540*x1*x6^2 + 5512*x1*x6 + 996*x1 + 1054*x2^3 + 7972*x2^2*x3 - 1870*x2^2*x5 - 8140*x2^2 + 18260*x2*x3*x5 - 3224*x2*x4^2 - 5998*x2*x4*x5 - 3080*x2*x4*x6 - 94*x2*x5^2 - 28412*x2*x5 - 1922*x2*x6 + 72*x3*x4^2 - 2016*x3*x4 + 13778*x3*x5^2 + 5408*x3*x6^2 + 14112*x3 + 912*x4^3 + 1140*x4^2*x5 - 4160*x4^2*x6 - 12768*x4^2 - 10790*x4*x5^2 - 4648*x4*x5*x6 - 16956*x4*x5 + 13944*x5 + 4680*x6^3 - 8632*x6,
    960*x1*x3*x4 + 3536*x1*x4*x5 - 8480*x1*x4 + 6372*x2^3 + 2296*x2^2*x4 + 6136*x2^2 - 6448*x2*x3*x4 - 5998*x2*x3*x5 - 3080*x2*x3*x6 + 29184*x2*x4*x5 + 11232*x2*x4 + 20450*x2*x5^2 + 952*x2*x5*x6 + 5624*x2*x5 + 2376*x2*x6^2 + 4144*x2*x6 - 3348*x2 + 72*x3^2*x4 - 1008*x3^2 + 2736*x3*x4^2 + 2280*x3*x4*x5 - 8320*x3*x4*x6 - 25536*x3*x4 - 10790*x3*x5^2 - 4648*x3*x5*x6 - 16956*x3*x5 + 40320*x4^3 + 43320*x4^2*x5 + 17348*x4*x5^2 + 7280*x4*x5*x6 - 25232*x4*x5 - 5632*x4*x6^2 + 6448*x4*x6 + 18688*x4 - 15770*x5^2 - 3848*x5 + 2288*x6^2 - 3224,
    578*x1^2*x5 - 578*x1*x2^2 - 1054*x1*x2*x3 + 1768*x1*x4^2 - 4488*x1*x5^2 + 1054*x1*x6 - 1870*x2^2*x3 + 22002*x2^2*x5 - 1850*x2^2 + 9130*x2*x3^2 - 5998*x2*x3*x4 - 188*x2*x3*x5 - 28412*x2*x3 + 14592*x2*x4^2 + 40900*x2*x4*x5 + 952*x2*x4*x6 + 5624*x2*x4 - 31872*x2*x5 + 13778*x3^2*x5 + 1140*x3*x4^2 - 21580*x3*x4*x5 - 4648*x3*x4*x6 - 16956*x3*x4 + 13944*x3 + 14440*x4^3 + 17348*x4^2*x5 + 3640*x4^2*x6 - 12616*x4^2 - 31540*x4*x5 - 3848*x4 + 7744*x5^3 - 5456*x5*x6 + 16516*x5 - 1628*x6^2 + 2294,
    -624*x1*x3^2 - 1080*x1*x3*x6 + 5512*x1*x3 + 1054*x1*x5 + 9540*x1*x6 + 5192*x2^2*x6 - 1054*x2^2 - 3080*x2*x3*x4 - 1922*x2*x3 + 952*x2*x4*x5 + 4752*x2*x4*x6 + 4144*x2*x4 + 5408*x3^2*x6 - 4160*x3*x4^2 - 4648*x3*x4*x5 + 14040*x3*x6^2 - 8632*x3 + 3640*x4^2*x5 - 5632*x4^2*x6 + 3224*x4^2 + 4576*x4*x6 - 2728*x5^2 - 3256*x5*x6 + 10036*x6^3 - 15746*x6,
    36*x1^2*x3^2 - 636*x1^2*x3 + 289*x1^2*x5^2 + 2809*x1^2 - 578*x1*x2^2*x5 - 1054*x1*x2*x3*x5 - 624*x1*x3^2*x6 + 480*x1*x3*x4^2 - 540*x1*x3*x6^2 + 5512*x1*x3*x6 + 996*x1*x3 + 1768*x1*x4^2*x5 - 4240*x1*x4^2 - 1496*x1*x5^3 + 1054*x1*x5*x6 + 4770*x1*x6^2 - 8798*x1 + 3770*x2^4 + 1054*x2^3*x3 + 6372*x2^3*x4 + 3986*x2^2*x3^2 - 1870*x2^2*x3*x5 - 8140*x2^2*x3 + 1148*x2^2*x4^2 + 6136*x2^2*x4 + 11001*x2^2*x5^2 - 1850*x2^2*x5 + 2596*x2^2*x6^2 - 1054*x2^2*x6 + 1818*x2^2 + 9130*x2*x3^2*x5 - 3224*x2*x3*x4^2 - 5998*x2*x3*x4*x5 - 3080*x2*x3*x4*x6 - 94*x2*x3*x5^2 - 28412*x2*x3*x5 - 1922*x2*x3*x6 + 14592*x2*x4^2*x5 + 5616*x2*x4^2 + 20450*x2*x4*x5^2 + 952*x2*x4*x5*x6 + 5624*x2*x4*x5 + 2376*x2*x4*x6^2 + 4144*x2*x4*x6 - 3348*x2*x4 - 15936*x2*x5^2 + 36*x3^2*x4^2 - 1008*x3^2*x4 + 6889*x3^2*x5^2 + 2704*x3^2*x6^2 + 7056*x3^2 + 912*x3*x4^3 + 1140*x3*x4^2*x5 - 4160*x3*x4^2*x6 - 12768*x3*x4^2 - 10790*x3*x4*x5^2 - 4648*x3*x4*x5*x6 - 16956*x3*x4*x5 + 13944*x3*x5 + 4680*x3*x6^3 - 8632*x3*x6 + 10080*x4^4 + 14440*x4^3*x5 + 8674*x4^2*x5^2 + 3640*x4^2*x5*x6 - 12616*x4^2*x5 - 2816*x4^2*x6^2 + 3224*x4^2*x6 + 9344*x4^2 - 15770*x4*x5^2 - 3848*x4*x5 + 2288*x4*x6^2 - 3224*x4 + 1936*x5^4 - 2728*x5^2*x6 + 8258*x5^2 - 1628*x5*x6^2 + 2294*x5 + 2509*x6^4 - 7873*x6^2 + 7850
]

# Reducing the rational system is far cheaper than building the same
# polynomials again over the prime field.
prime = Int32(AlgebraicSolving.Nemo.rand_bits_prime(ZZ, 31))
Rp, _ = polynomial_ring(GF(prime), ["x1", "x2", "x3", "x4", "x5", "x6"])
Fp = [AlgebraicSolving.reduce_mod_p(f, Rp) for f in F]

# Decompose a tiny system first so that the timings below measure the
# computation rather than compilation. The solver specialises on the number of
# variables, so the warm up uses as many of them as the benchmark itself.
function warmup(K)
    _, w = polynomial_ring(K, ["w$(i)" for i in 1:6])
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
        println(io, join((NAME, 6, 6,
                          round(time_qq, digits = 6),
                          round(time_gf, digits = 6),
                          round(ratio, digits = 4)), ","))
    end
end
println("appended to $(csv)")
