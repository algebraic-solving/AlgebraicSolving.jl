# Benchmark for `equidimensional_decomposition` on the "soscritpoints_deg_2_nvars_6_s_4" system.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AlgebraicSolving
using FileWatching.Pidfile: mkpidlock

const NAME = "soscritpoints_deg_2_nvars_6_s_4"

R, (x1, x2, x3, x4, x5, x6) = polynomial_ring(QQ, ["x1", "x2", "x3", "x4", "x5", "x6"])

F = [
    4524*x1*x2*x5 + 2496*x1*x2 - 344*x1*x3*x5 - 228*x1*x3*x6 + 3608*x1*x4*x5 - 7872*x1*x5^2 + 960*x1*x5 + 2040*x1 + 676*x2^3 + 4524*x2*x3*x4 + 3276*x2*x3*x6 + 28240*x2*x5^2 + 19608*x2*x5*x6 + 6498*x2*x6^2 + 450*x2 - 510*x3^2 - 6020*x3*x5^2 + 5814*x3*x5*x6 + 6498*x3*x6^2 + 10492*x4*x5 + 4344*x4*x6 - 7544*x5^2*x6 + 8036*x5*x6^2 + 8200*x5*x6 + 1800*x6^2,
    8*x1^2*x3 - 344*x1*x2*x5 - 228*x1*x2*x6 - 1896*x1*x3*x5 - 456*x1*x3*x6 - 4624*x1*x3 + 15138*x1*x4*x5 + 8108*x1*x4 + 10962*x1*x5*x6 + 6048*x1*x6 + 2262*x2^2*x4 + 1638*x2^2*x6 - 1020*x2*x3 - 6020*x2*x5^2 + 5814*x2*x5*x6 + 6498*x2*x6^2 + 1156*x3^3 + 15138*x3*x4^2 + 27840*x3*x4*x6 + 2450*x3*x5^2 - 7980*x3*x5*x6 + 10356*x3*x6^2 - 4270*x4*x5 + 6954*x4*x6,
    968*x1^2*x4 - 2112*x1^2*x5 + 3608*x1*x2*x5 + 15138*x1*x3*x5 + 8108*x1*x3 - 7592*x1*x5*x6 + 2156*x1*x6^2 - 9632*x1*x6 + 2262*x2^2*x3 + 10492*x2*x5 + 4344*x2*x6 + 15138*x3^2*x4 + 13920*x3^2*x6 - 4270*x3*x5 + 6954*x3*x6 + 15138*x4*x6^2 + 7442*x4 - 10440*x6^3,
    -2112*x1^2*x4 + 21794*x1^2*x5 + 12704*x1^2 + 2262*x1*x2^2 - 344*x1*x2*x3 + 3608*x1*x2*x4 - 15744*x1*x2*x5 + 960*x1*x2 - 948*x1*x3^2 + 15138*x1*x3*x4 + 10962*x1*x3*x6 - 7592*x1*x4*x6 + 8832*x1*x5*x6 - 864*x1*x6^2 - 4800*x1*x6 + 28240*x2^2*x5 + 9804*x2^2*x6 - 12040*x2*x3*x5 + 5814*x2*x3*x6 + 10492*x2*x4 - 15088*x2*x5*x6 + 8036*x2*x6^2 + 8200*x2*x6 + 2450*x3^2*x5 - 3990*x3^2*x6 - 4270*x3*x4 + 4232*x5*x6^2 - 4508*x6^3 - 4600*x6^2,
    -228*x1*x2*x3 - 228*x1*x3^2 + 10962*x1*x3*x5 + 6048*x1*x3 - 7592*x1*x4*x5 + 4312*x1*x4*x6 - 9632*x1*x4 + 4416*x1*x5^2 - 1728*x1*x5*x6 - 4800*x1*x5 + 16320*x1*x6 + 1638*x2^2*x3 + 9804*x2^2*x5 + 6498*x2^2*x6 + 5814*x2*x3*x5 + 12996*x2*x3*x6 + 4344*x2*x4 - 7544*x2*x5^2 + 16072*x2*x5*x6 + 8200*x2*x5 + 3600*x2*x6 + 13920*x3^2*x4 - 3990*x3^2*x5 + 10356*x3^2*x6 + 6954*x3*x4 + 15138*x4^2*x6 - 31320*x4*x6^2 + 4232*x5^2*x6 - 13524*x5*x6^2 - 9200*x5*x6 + 24004*x6^3 + 14700*x6^2 + 5000*x6,
    4*x1^2*x3^2 + 484*x1^2*x4^2 - 2112*x1^2*x4*x5 + 10897*x1^2*x5^2 + 12704*x1^2*x5 + 6928*x1^2 + 2262*x1*x2^2*x5 + 1248*x1*x2^2 - 344*x1*x2*x3*x5 - 228*x1*x2*x3*x6 + 3608*x1*x2*x4*x5 - 7872*x1*x2*x5^2 + 960*x1*x2*x5 + 2040*x1*x2 - 948*x1*x3^2*x5 - 228*x1*x3^2*x6 - 2312*x1*x3^2 + 15138*x1*x3*x4*x5 + 8108*x1*x3*x4 + 10962*x1*x3*x5*x6 + 6048*x1*x3*x6 - 7592*x1*x4*x5*x6 + 2156*x1*x4*x6^2 - 9632*x1*x4*x6 + 4416*x1*x5^2*x6 - 864*x1*x5*x6^2 - 4800*x1*x5*x6 + 8160*x1*x6^2 + 169*x2^4 + 2262*x2^2*x3*x4 + 1638*x2^2*x3*x6 + 14120*x2^2*x5^2 + 9804*x2^2*x5*x6 + 3249*x2^2*x6^2 + 225*x2^2 - 510*x2*x3^2 - 6020*x2*x3*x5^2 + 5814*x2*x3*x5*x6 + 6498*x2*x3*x6^2 + 10492*x2*x4*x5 + 4344*x2*x4*x6 - 7544*x2*x5^2*x6 + 8036*x2*x5*x6^2 + 8200*x2*x5*x6 + 1800*x2*x6^2 + 289*x3^4 + 7569*x3^2*x4^2 + 13920*x3^2*x4*x6 + 1225*x3^2*x5^2 - 3990*x3^2*x5*x6 + 5178*x3^2*x6^2 - 4270*x3*x4*x5 + 6954*x3*x4*x6 + 7569*x4^2*x6^2 + 3721*x4^2 - 10440*x4*x6^3 + 2116*x5^2*x6^2 - 4508*x5*x6^3 - 4600*x5*x6^2 + 6001*x6^4 + 4900*x6^3 + 2500*x6^2
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
