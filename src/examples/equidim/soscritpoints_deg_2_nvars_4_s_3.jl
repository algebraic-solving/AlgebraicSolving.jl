# Benchmark for `equidimensional_decomposition` on the "soscritpoints_deg_2_nvars_4_s_3" system.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AlgebraicSolving
using FileWatching.Pidfile: mkpidlock

const NAME = "soscritpoints_deg_2_nvars_4_s_3"

R, (x1, x2, x3, x4) = polynomial_ring(QQ, ["x1", "x2", "x3", "x4"])

F = [
    -5510*x1^3 + 18050*x1^2*x2 - 638*x1^2*x3 + 4184*x1*x2*x3 - 27600*x1*x2 - 9310*x1*x3^2 - 8820*x1*x3*x4 + 1090*x1*x3 + 11100*x1*x4 - 1520*x1 + 33860*x2^3 - 3312*x2^2*x3 - 40518*x2^2*x4 + 570*x2^2 - 26182*x2*x3^2 - 11840*x2*x3*x4 + 16890*x2*x4^2 + 20900*x2*x4 + 18358*x2 - 214*x3^3 + 10066*x3^2*x4 + 5476*x3*x4^2 - 176*x3 - 3080*x4^3 - 5320*x4^2 + 8470*x4 + 14630,
    -638*x1^2*x2 + 5686*x1^2*x3 + 2726*x1^2*x4 + 2092*x1*x2^2 - 18620*x1*x2*x3 - 8820*x1*x2*x4 + 1090*x1*x2 + 21600*x1*x3 - 56*x1*x4^2 + 5550*x1*x4 + 154*x1 - 1104*x2^3 - 26182*x2^2*x3 - 5920*x2^2*x4 - 642*x2*x3^2 + 20132*x2*x3*x4 + 5476*x2*x4^2 - 176*x2 + 30340*x3^3 + 29802*x3^2*x4 + 7156*x3*x4^2 + 1568*x3 + 752*x4,
    2726*x1^2*x3 - 8820*x1*x2*x3 + 11100*x1*x2 - 112*x1*x3*x4 + 5550*x1*x3 - 13506*x2^3 - 5920*x2^2*x3 + 16890*x2^2*x4 + 10450*x2^2 + 10066*x2*x3^2 + 10952*x2*x3*x4 - 9240*x2*x4^2 - 10640*x2*x4 + 8470*x2 + 9934*x3^3 + 7156*x3^2*x4 + 752*x3 + 3136*x4^3 - 8624*x4,
    841*x1^4 - 5510*x1^3*x2 + 9025*x1^2*x2^2 - 638*x1^2*x2*x3 + 2843*x1^2*x3^2 + 2726*x1^2*x3*x4 + 6089*x1^2 + 2092*x1*x2^2*x3 - 13800*x1*x2^2 - 9310*x1*x2*x3^2 - 8820*x1*x2*x3*x4 + 1090*x1*x2*x3 + 11100*x1*x2*x4 - 1520*x1*x2 + 10800*x1*x3^2 - 56*x1*x3*x4^2 + 5550*x1*x3*x4 + 154*x1*x3 + 8465*x2^4 - 1104*x2^3*x3 - 13506*x2^3*x4 + 190*x2^3 - 13091*x2^2*x3^2 - 5920*x2^2*x3*x4 + 8445*x2^2*x4^2 + 10450*x2^2*x4 + 9179*x2^2 - 214*x2*x3^3 + 10066*x2*x3^2*x4 + 5476*x2*x3*x4^2 - 176*x2*x3 - 3080*x2*x4^3 - 5320*x2*x4^2 + 8470*x2*x4 + 14630*x2 + 7585*x3^4 + 9934*x3^3*x4 + 3578*x3^2*x4^2 + 784*x3^2 + 752*x3*x4 + 784*x4^4 - 4312*x4^2 + 5993
]

# Reducing the rational system is far cheaper than building the same
# polynomials again over the prime field.
prime = Int32(AlgebraicSolving.Nemo.rand_bits_prime(ZZ, 31))
Rp, _ = polynomial_ring(GF(prime), ["x1", "x2", "x3", "x4"])
Fp = [AlgebraicSolving.reduce_mod_p(f, Rp) for f in F]

# Decompose a tiny system first so that the timings below measure the
# computation rather than compilation. The solver specialises on the number of
# variables, so the warm up uses as many of them as the benchmark itself.
function warmup(K)
    _, w = polynomial_ring(K, ["w$(i)" for i in 1:4])
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
        println(io, join((NAME, 4, 4,
                          round(time_qq, digits = 6),
                          round(time_gf, digits = 6),
                          round(ratio, digits = 4)), ","))
    end
end
println("appended to $(csv)")
