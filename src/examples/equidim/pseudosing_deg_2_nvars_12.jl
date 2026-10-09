# Benchmark for `equidimensional_decomposition` on the "pseudosing_deg_2_nvars_12" system.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AlgebraicSolving
using FileWatching.Pidfile: mkpidlock

const NAME = "pseudosing_deg_2_nvars_12"

R, (x1, x2, x3, x4, x5, x6, y1, y2, y3, y4, y5, y6) = polynomial_ring(QQ, ["x1", "x2", "x3", "x4", "x5", "x6", "y1", "y2", "y3", "y4", "y5", "y6"])

F = [
    x6 - y6,
    x5 - y5,
    -46*x1^2 - 33*x1*x2 + 87*x1*x3 - 34*x1*x4 + 40*x1*x5 + 77*x1*x6 + x1 - 10*x2^2 - 65*x2*x3 - 85*x2*x4 + 54*x2*x5 + 18*x2 + 52*x3^2 + 36*x3*x4 + 91*x3*x5 - 22*x3*x6 + 51*x3 - 27*x4^2 + 50*x4*x5 + 60*x4*x6 - 91*x4 - 47*x5^2 - 97*x5*x6 - 2*x5 - 31*x6^2 + 25*x6 + 31,
    -27*x1^2 + 65*x1*x2 + 88*x1*x3 + 10*x1*x4 - 6*x1*x5 + 80*x1*x6 - 84*x1 + 57*x2^2 - 49*x2*x3 + 31*x2*x4 + 73*x2*x5 + 95*x2*x6 + 68*x2 - 29*x3^2 + 5*x3*x4 - 26*x3*x5 - 51*x3*x6 + 88*x3 + 97*x4^2 - 67*x4*x5 + 58*x4*x6 + 29*x4 + 37*x5^2 + 5*x5*x6 - 36*x5 - 57*x6^2 + 85*x6 + 80,
    90*x1^2 + 74*x1*x2 + 27*x1*x3 + 9*x1*x4 - 91*x1*x5 + 81*x1*x6 + 65*x1 - 12*x2^2 + 78*x2*x3 + 5*x2*x4 - 63*x2*x5 - 5*x2*x6 + 36*x2 - 8*x3^2 + 30*x3*x4 - 3*x3*x5 - 56*x3*x6 - 91*x3 - 70*x4^2 + 42*x4*x5 + 9*x4*x6 - 21*x4 - 27*x5^2 - 79*x5*x6 - 22*x5 - 51*x6^2 + 16*x6 - 85,
    -44*x1^2 - 31*x1*x2 + 45*x1*x3 + 49*x1*x4 - 58*x1*x5 + 49*x1*x6 + x1 - 95*x2^2 + 86*x2*x3 - 97*x2*x4 - 14*x2*x5 + 83*x2*x6 - 96*x2 - 8*x3^2 - 54*x3*x4 + 62*x3*x5 + 96*x3*x6 - 51*x3 + 89*x4^2 + 14*x4*x5 - 79*x4*x6 - 58*x4 - 95*x5^2 + 61*x5*x6 - 2*x5 + 86*x6^2 + 57*x6 - 35,
    57*x1^2 + 28*x1*x2 + 63*x1*x3 + 21*x1*x4 - 71*x1*x5 - 66*x1*x6 - 34*x1 + 72*x2^2 - 40*x2*x3 - 68*x2*x4 - 15*x2*x5 - 32*x2*x6 + 17*x2 + 87*x3^2 - 60*x3*x4 + 7*x3*x5 - 61*x3*x6 + 45*x3 - 50*x4^2 - 22*x4*x5 + 48*x4*x6 - 82*x4 + 46*x5^2 - 49*x5*x6 - 66*x5 + 18*x6^2 + 16*x6 - 22,
    -46*y1^2 - 33*y1*y2 + 87*y1*y3 - 34*y1*y4 + 40*y1*y5 + 77*y1*y6 + y1 - 10*y2^2 - 65*y2*y3 - 85*y2*y4 + 54*y2*y5 + 18*y2 + 52*y3^2 + 36*y3*y4 + 91*y3*y5 - 22*y3*y6 + 51*y3 - 27*y4^2 + 50*y4*y5 + 60*y4*y6 - 91*y4 - 47*y5^2 - 97*y5*y6 - 2*y5 - 31*y6^2 + 25*y6 + 31,
    -27*y1^2 + 65*y1*y2 + 88*y1*y3 + 10*y1*y4 - 6*y1*y5 + 80*y1*y6 - 84*y1 + 57*y2^2 - 49*y2*y3 + 31*y2*y4 + 73*y2*y5 + 95*y2*y6 + 68*y2 - 29*y3^2 + 5*y3*y4 - 26*y3*y5 - 51*y3*y6 + 88*y3 + 97*y4^2 - 67*y4*y5 + 58*y4*y6 + 29*y4 + 37*y5^2 + 5*y5*y6 - 36*y5 - 57*y6^2 + 85*y6 + 80,
    90*y1^2 + 74*y1*y2 + 27*y1*y3 + 9*y1*y4 - 91*y1*y5 + 81*y1*y6 + 65*y1 - 12*y2^2 + 78*y2*y3 + 5*y2*y4 - 63*y2*y5 - 5*y2*y6 + 36*y2 - 8*y3^2 + 30*y3*y4 - 3*y3*y5 - 56*y3*y6 - 91*y3 - 70*y4^2 + 42*y4*y5 + 9*y4*y6 - 21*y4 - 27*y5^2 - 79*y5*y6 - 22*y5 - 51*y6^2 + 16*y6 - 85,
    -44*y1^2 - 31*y1*y2 + 45*y1*y3 + 49*y1*y4 - 58*y1*y5 + 49*y1*y6 + y1 - 95*y2^2 + 86*y2*y3 - 97*y2*y4 - 14*y2*y5 + 83*y2*y6 - 96*y2 - 8*y3^2 - 54*y3*y4 + 62*y3*y5 + 96*y3*y6 - 51*y3 + 89*y4^2 + 14*y4*y5 - 79*y4*y6 - 58*y4 - 95*y5^2 + 61*y5*y6 - 2*y5 + 86*y6^2 + 57*y6 - 35,
    57*y1^2 + 28*y1*y2 + 63*y1*y3 + 21*y1*y4 - 71*y1*y5 - 66*y1*y6 - 34*y1 + 72*y2^2 - 40*y2*y3 - 68*y2*y4 - 15*y2*y5 - 32*y2*y6 + 17*y2 + 87*y3^2 - 60*y3*y4 + 7*y3*y5 - 61*y3*y6 + 45*y3 - 50*y4^2 - 22*y4*y5 + 48*y4*y6 - 82*y4 + 46*y5^2 - 49*y5*y6 - 66*y5 + 18*y6^2 + 16*y6 - 22
]

# Reducing the rational system is far cheaper than building the same
# polynomials again over the prime field.
prime = Int32(AlgebraicSolving.Nemo.rand_bits_prime(ZZ, 31))
Rp, _ = polynomial_ring(GF(prime), ["x1", "x2", "x3", "x4", "x5", "x6", "y1", "y2", "y3", "y4", "y5", "y6"])
Fp = [AlgebraicSolving.reduce_mod_p(f, Rp) for f in F]

# Decompose a tiny system first so that the timings below measure the
# computation rather than compilation. The solver specialises on the number of
# variables, so the warm up uses as many of them as the benchmark itself.
function warmup(K)
    _, w = polynomial_ring(K, ["w$(i)" for i in 1:12])
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
        println(io, join((NAME, 12, 12,
                          round(time_qq, digits = 6),
                          round(time_gf, digits = 6),
                          round(ratio, digits = 4)), ","))
    end
end
println("appended to $(csv)")
