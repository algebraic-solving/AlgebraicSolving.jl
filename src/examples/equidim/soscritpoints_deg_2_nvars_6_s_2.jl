# Benchmark for `equidimensional_decomposition` on the "soscritpoints_deg_2_nvars_6_s_2" system.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AlgebraicSolving
using FileWatching.Pidfile: mkpidlock

const NAME = "soscritpoints_deg_2_nvars_6_s_2"

R, (x1, x2, x3, x4, x5, x6) = polynomial_ring(QQ, ["x1", "x2", "x3", "x4", "x5", "x6"])

F = [
    1682*x1^2*x2 + 2146*x1^2*x3 + 4368*x1*x2*x4 - 3308*x1*x2*x5 - 3120*x1*x2*x6 + 370*x1*x3*x5 - 3886*x1*x4 - 2088*x1*x6^2 + 3364*x1*x6 + 576*x2^3 - 3744*x2*x3*x4 - 240*x2*x4*x6 + 50*x2*x5^2 - 670*x4*x5 - 360*x5*x6^2 + 580*x5*x6,
    2146*x1^2*x2 + 2738*x1^2*x3 + 370*x1*x2*x5 - 14196*x1*x4^2 + 12636*x1*x4*x5 + 10140*x1*x4*x6 - 4958*x1*x4 - 2664*x1*x6^2 + 4292*x1*x6 - 1872*x2^2*x4 + 12168*x3*x4^2 + 780*x4^2*x6,
    16562*x1^2*x4 - 14742*x1^2*x5 - 11830*x1^2*x6 + 2184*x1*x2^2 - 3886*x1*x2 - 28392*x1*x3*x4 + 12636*x1*x3*x5 + 10140*x1*x3*x6 - 4958*x1*x3 - 1820*x1*x4*x6 + 810*x1*x5*x6 + 650*x1*x6^2 - 1872*x2^2*x3 - 120*x2^2*x6 - 670*x2*x5 + 12168*x3^2*x4 + 1560*x3*x4*x6 + 50*x4*x6^2 + 8978*x4 + 4824*x6^2 - 7772*x6,
    -14742*x1^2*x4 + 13122*x1^2*x5 + 10530*x1^2*x6 - 1654*x1*x2^2 + 370*x1*x2*x3 + 12636*x1*x3*x4 + 810*x1*x4*x6 + 50*x2^2*x5 - 670*x2*x4 - 360*x2*x6^2 + 580*x2*x6,
    -11830*x1^2*x4 + 10530*x1^2*x5 + 8450*x1^2*x6 - 1560*x1*x2^2 - 4176*x1*x2*x6 + 3364*x1*x2 + 10140*x1*x3*x4 - 5328*x1*x3*x6 + 4292*x1*x3 - 910*x1*x4^2 + 810*x1*x4*x5 + 1300*x1*x4*x6 - 120*x2^2*x4 - 720*x2*x5*x6 + 580*x2*x5 + 780*x3*x4^2 + 50*x4^2*x6 + 9648*x4*x6 - 7772*x4 + 5184*x6^3 - 12528*x6^2 + 6728*x6,
    841*x1^2*x2^2 + 2146*x1^2*x2*x3 + 1369*x1^2*x3^2 + 8281*x1^2*x4^2 - 14742*x1^2*x4*x5 - 11830*x1^2*x4*x6 + 6561*x1^2*x5^2 + 10530*x1^2*x5*x6 + 4225*x1^2*x6^2 + 2184*x1*x2^2*x4 - 1654*x1*x2^2*x5 - 1560*x1*x2^2*x6 + 370*x1*x2*x3*x5 - 3886*x1*x2*x4 - 2088*x1*x2*x6^2 + 3364*x1*x2*x6 - 14196*x1*x3*x4^2 + 12636*x1*x3*x4*x5 + 10140*x1*x3*x4*x6 - 4958*x1*x3*x4 - 2664*x1*x3*x6^2 + 4292*x1*x3*x6 - 910*x1*x4^2*x6 + 810*x1*x4*x5*x6 + 650*x1*x4*x6^2 + 144*x2^4 - 1872*x2^2*x3*x4 - 120*x2^2*x4*x6 + 25*x2^2*x5^2 - 670*x2*x4*x5 - 360*x2*x5*x6^2 + 580*x2*x5*x6 + 6084*x3^2*x4^2 + 780*x3*x4^2*x6 + 25*x4^2*x6^2 + 4489*x4^2 + 4824*x4*x6^2 - 7772*x4*x6 + 1296*x6^4 - 4176*x6^3 + 3364*x6^2
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
