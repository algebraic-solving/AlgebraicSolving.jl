# Benchmark for `equidimensional_decomposition` on the "soscritpoints_deg_2_nvars_6_s_3" system.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AlgebraicSolving
using FileWatching.Pidfile: mkpidlock

const NAME = "soscritpoints_deg_2_nvars_6_s_3"

R, (x1, x2, x3, x4, x5, x6) = polynomial_ring(QQ, ["x1", "x2", "x3", "x4", "x5", "x6"])

F = [
    512*x1^2*x2 - 8160*x1*x2^2 - 2816*x1*x2*x5 - 992*x1*x3*x5 - 1328*x1*x3 + 1440*x1*x4^2 + 1568*x1*x6^2 + 28900*x2^3 + 22440*x2^2*x5 + 10540*x2*x3*x5 - 15300*x2*x4^2 + 3872*x2*x5^2 - 16660*x2*x6^2 + 13778*x2 + 2728*x3*x5^2 - 8964*x3*x5 + 10292*x3*x6 - 3960*x4^2*x5 - 15936*x4 + 15936*x5^2 - 4312*x5*x6^2,
    128*x1^2*x3 - 992*x1*x2*x5 - 1328*x1*x2 + 1728*x1*x3*x5 - 1984*x1*x3*x6 + 1536*x1*x4 - 1536*x1*x5^2 + 5270*x2^2*x5 + 2728*x2*x5^2 - 8964*x2*x5 + 10292*x2*x6 + 7754*x3*x5^2 - 13392*x3*x5*x6 + 7688*x3*x6^2 - 2790*x4^2*x5 + 10368*x4*x5 - 11904*x4*x6 - 10368*x5^3 + 11904*x5^2*x6 - 3038*x5*x6^2,
    -5880*x1^2 + 2880*x1*x2*x4 + 1536*x1*x3 - 1260*x1*x5 + 12740*x1 - 15300*x2^2*x4 - 7920*x2*x4*x5 - 15936*x2 - 5580*x3*x4*x5 + 10368*x3*x5 - 11904*x3*x6 + 8100*x4^3 + 8820*x4*x6^2 + 28232*x4 - 18432*x5^2 + 2940*x6^2 + 7840,
    756*x1^3 + 162*x1^2*x5 - 1638*x1^2 - 1408*x1*x2^2 - 992*x1*x2*x3 + 864*x1*x3^2 - 3072*x1*x3*x5 - 1260*x1*x4 - 378*x1*x6^2 - 1008*x1 + 7480*x2^3 + 5270*x2^2*x3 + 3872*x2^2*x5 + 5456*x2*x3*x5 - 8964*x2*x3 - 3960*x2*x4^2 + 31872*x2*x5 - 4312*x2*x6^2 + 7754*x3^2*x5 - 6696*x3^2*x6 - 2790*x3*x4^2 + 10368*x3*x4 - 31104*x3*x5^2 + 23808*x3*x5*x6 - 3038*x3*x6^2 - 36864*x4*x5 + 36864*x5^3,
    -3528*x1^2*x6 + 3136*x1*x2*x6 - 992*x1*x3^2 - 756*x1*x5*x6 + 7644*x1*x6 - 16660*x2^2*x6 + 10292*x2*x3 - 8624*x2*x5*x6 - 6696*x3^2*x5 + 7688*x3^2*x6 - 11904*x3*x4 + 11904*x3*x5^2 - 6076*x3*x5*x6 + 8820*x4^2*x6 + 5880*x4*x6 + 11368*x6^3 + 4704*x6,
    1764*x1^4 + 756*x1^3*x5 - 7644*x1^3 + 256*x1^2*x2^2 + 64*x1^2*x3^2 - 5880*x1^2*x4 + 81*x1^2*x5^2 - 1638*x1^2*x5 - 1764*x1^2*x6^2 + 3577*x1^2 - 2720*x1*x2^3 - 1408*x1*x2^2*x5 - 992*x1*x2*x3*x5 - 1328*x1*x2*x3 + 1440*x1*x2*x4^2 + 1568*x1*x2*x6^2 + 864*x1*x3^2*x5 - 992*x1*x3^2*x6 + 1536*x1*x3*x4 - 1536*x1*x3*x5^2 - 1260*x1*x4*x5 + 12740*x1*x4 - 378*x1*x5*x6^2 - 1008*x1*x5 + 3822*x1*x6^2 + 10192*x1 + 7225*x2^4 + 7480*x2^3*x5 + 5270*x2^2*x3*x5 - 7650*x2^2*x4^2 + 1936*x2^2*x5^2 - 8330*x2^2*x6^2 + 6889*x2^2 + 2728*x2*x3*x5^2 - 8964*x2*x3*x5 + 10292*x2*x3*x6 - 3960*x2*x4^2*x5 - 15936*x2*x4 + 15936*x2*x5^2 - 4312*x2*x5*x6^2 + 3877*x3^2*x5^2 - 6696*x3^2*x5*x6 + 3844*x3^2*x6^2 - 2790*x3*x4^2*x5 + 10368*x3*x4*x5 - 11904*x3*x4*x6 - 10368*x3*x5^3 + 11904*x3*x5^2*x6 - 3038*x3*x5*x6^2 + 2025*x4^4 + 4410*x4^2*x6^2 + 14116*x4^2 - 18432*x4*x5^2 + 2940*x4*x6^2 + 7840*x4 + 9216*x5^4 + 2842*x6^4 + 2352*x6^2 + 3136
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
