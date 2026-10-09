# Benchmark for `equidimensional_decomposition` on the "soscritpoints_deg_2_nvars_5_s_4" system.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AlgebraicSolving
using FileWatching.Pidfile: mkpidlock

const NAME = "soscritpoints_deg_2_nvars_5_s_4"

R, (x1, x2, x3, x4, x5) = polynomial_ring(QQ, ["x1", "x2", "x3", "x4", "x5"])

F = [
    18130*x1^2*x2 + 6202*x1^2*x3 - 1848*x1^2*x5 + 35308*x1*x2*x4 + 728*x1*x2*x5 - 18200*x1*x2 + 1176*x1*x3^2 + 13644*x1*x3*x4 + 188*x1*x3*x5 - 12086*x1*x3 - 1960*x1*x4*x5 + 4480*x1*x5 - 5916*x1 + 29476*x2*x4^2 + 776*x2*x4*x5 - 2756*x2*x4 + 8*x2*x5^2 - 400*x2*x5 + 13810*x2 + 13870*x3^2*x4 + 10830*x3^2 - 11708*x3*x4 - 240*x3*x5 + 3280*x3 + 9928*x4^2*x5 - 7154*x4^2 + 8432*x4*x5 - 10822*x4,
    6202*x1^2*x2 + 9868*x1^2*x3 + 2772*x1^2*x5 + 2352*x1*x2*x3 + 13644*x1*x2*x4 + 188*x1*x2*x5 - 12086*x1*x2 + 12378*x1*x3^2 - 2772*x1*x3*x5 - 11280*x1*x3 + 7156*x1*x4*x5 - 2864*x1*x4 - 6720*x1*x5 + 6960*x1 + 27740*x2*x3*x4 + 21660*x2*x3 - 11708*x2*x4 - 240*x2*x5 + 3280*x2 + 37864*x3^3 + 2*x3*x4^2 + 22900*x3*x4*x5 - 18460*x3*x4 + 6720*x3*x5 + 10400*x3 - 20*x4^2*x5 + 154*x4^2 - 800*x4*x5 + 6160*x4,
    17654*x1*x2^2 + 13644*x1*x2*x3 - 1960*x1*x2*x5 + 7156*x1*x3*x5 - 2864*x1*x3 + 2310*x1*x5^2 - 1740*x1*x5 + 13398*x1 + 29476*x2^2*x4 + 388*x2^2*x5 - 1378*x2^2 + 13870*x2*x3^2 - 11708*x2*x3 + 19856*x2*x4*x5 - 14308*x2*x4 + 8432*x2*x5 - 10822*x2 + 2*x3^2*x4 + 11450*x3^2*x5 - 9230*x3^2 - 40*x3*x4*x5 + 308*x3*x4 - 800*x3*x5 + 6160*x3 + 11898*x4*x5^2 - 16408*x4*x5 + 16660*x4 - 5600*x5^2,
    -1848*x1^2*x2 + 2772*x1^2*x3 + 2178*x1^2*x5 + 364*x1*x2^2 + 188*x1*x2*x3 - 1960*x1*x2*x4 + 4480*x1*x2 - 1386*x1*x3^2 + 7156*x1*x3*x4 - 6720*x1*x3 + 4620*x1*x4*x5 - 1740*x1*x4 - 10560*x1*x5 + 388*x2^2*x4 + 8*x2^2*x5 - 200*x2^2 - 240*x2*x3 + 9928*x2*x4^2 + 8432*x2*x4 + 11450*x3^2*x4 + 3360*x3^2 - 20*x3*x4^2 - 800*x3*x4 + 11898*x4^2*x5 - 8204*x4^2 - 11200*x4*x5 + 12800*x5,
    9065*x1^2*x2^2 + 6202*x1^2*x2*x3 - 1848*x1^2*x2*x5 + 4934*x1^2*x3^2 + 2772*x1^2*x3*x5 + 1089*x1^2*x5^2 + 7569*x1^2 + 17654*x1*x2^2*x4 + 364*x1*x2^2*x5 - 9100*x1*x2^2 + 1176*x1*x2*x3^2 + 13644*x1*x2*x3*x4 + 188*x1*x2*x3*x5 - 12086*x1*x2*x3 - 1960*x1*x2*x4*x5 + 4480*x1*x2*x5 - 5916*x1*x2 + 4126*x1*x3^3 - 1386*x1*x3^2*x5 - 5640*x1*x3^2 + 7156*x1*x3*x4*x5 - 2864*x1*x3*x4 - 6720*x1*x3*x5 + 6960*x1*x3 + 2310*x1*x4*x5^2 - 1740*x1*x4*x5 + 13398*x1*x4 - 5280*x1*x5^2 + 14738*x2^2*x4^2 + 388*x2^2*x4*x5 - 1378*x2^2*x4 + 4*x2^2*x5^2 - 200*x2^2*x5 + 6905*x2^2 + 13870*x2*x3^2*x4 + 10830*x2*x3^2 - 11708*x2*x3*x4 - 240*x2*x3*x5 + 3280*x2*x3 + 9928*x2*x4^2*x5 - 7154*x2*x4^2 + 8432*x2*x4*x5 - 10822*x2*x4 + 9466*x3^4 + x3^2*x4^2 + 11450*x3^2*x4*x5 - 9230*x3^2*x4 + 3360*x3^2*x5 + 5200*x3^2 - 20*x3*x4^2*x5 + 154*x3*x4^2 - 800*x3*x4*x5 + 6160*x3*x4 + 5949*x4^2*x5^2 - 8204*x4^2*x5 + 8330*x4^2 - 5600*x4*x5^2 + 6400*x5^2
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
