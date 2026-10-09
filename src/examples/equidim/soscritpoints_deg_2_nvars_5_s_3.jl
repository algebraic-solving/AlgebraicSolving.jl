# Benchmark for `equidimensional_decomposition` on the "soscritpoints_deg_2_nvars_5_s_3" system.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AlgebraicSolving
using FileWatching.Pidfile: mkpidlock

const NAME = "soscritpoints_deg_2_nvars_5_s_3"

R, (x1, x2, x3, x4, x5) = polynomial_ring(QQ, ["x1", "x2", "x3", "x4", "x5"])

F = [
    21474*x1^2*x2 - 66*x1^2*x4 + 4760*x1^2 + 5616*x1*x3^2 - 3264*x1*x3*x4 - 5538*x1*x3 + 9840*x1*x4*x5 + 17472*x1*x5^2 + 7626*x1*x5 + 3528*x2*x4^2 + 4200*x2*x4 + 1250*x2 - 588*x3^2*x4 - 350*x3^2 + 8064*x4^2 + 4800*x4,
    1400*x1^2*x3 + 11232*x1*x2*x3 - 3264*x1*x2*x4 - 5538*x1*x2 + 15264*x1*x3*x4 - 1680*x1*x3*x5 + 576*x1*x3 - 8784*x1*x4 - 284*x1 - 1176*x2*x3*x4 - 700*x2*x3 + 20932*x3^3 - 30672*x3^2 + 578*x3*x4^2 - 2688*x3*x4 + 13536*x3*x5 + 10082*x3 - 850*x4^2*x5 - 3094*x4*x5^2 - 170*x4*x5 - 6674*x5,
    -66*x1^2*x2 + 5618*x1^2*x4 - 9388*x1^2 - 3264*x1*x2*x3 + 9840*x1*x2*x5 + 7632*x1*x3^2 - 8784*x1*x3 + 18352*x1*x5 + 3528*x2^2*x4 + 2100*x2^2 - 588*x2*x3^2 + 16128*x2*x4 + 4800*x2 + 578*x3^2*x4 - 1344*x3^2 - 1700*x3*x4*x5 - 3094*x3*x5^2 - 170*x3*x5 + 1250*x4*x5^2 + 18432*x4 + 4550*x5^3 + 250*x5^2,
    -6000*x1^3 + 7200*x1^2*x5 + 9840*x1*x2*x4 + 34944*x1*x2*x5 + 7626*x1*x2 - 840*x1*x3^2 + 18352*x1*x4 + 13468*x1*x5 + 558*x1 + 6768*x3^2 - 850*x3*x4^2 - 6188*x3*x4*x5 - 170*x3*x4 - 6674*x3 + 1250*x4^2*x5 + 13650*x4*x5^2 + 500*x4*x5 + 33124*x5^3 + 2730*x5^2 + 4468*x5,
    2500*x1^4 - 6000*x1^3*x5 + 10737*x1^2*x2^2 - 66*x1^2*x2*x4 + 4760*x1^2*x2 + 700*x1^2*x3^2 + 2809*x1^2*x4^2 - 9388*x1^2*x4 + 3600*x1^2*x5^2 + 1373*x1^2 + 5616*x1*x2*x3^2 - 3264*x1*x2*x3*x4 - 5538*x1*x2*x3 + 9840*x1*x2*x4*x5 + 17472*x1*x2*x5^2 + 7626*x1*x2*x5 + 7632*x1*x3^2*x4 - 840*x1*x3^2*x5 + 288*x1*x3^2 - 8784*x1*x3*x4 - 284*x1*x3 + 18352*x1*x4*x5 + 6734*x1*x5^2 + 558*x1*x5 + 1764*x2^2*x4^2 + 2100*x2^2*x4 + 625*x2^2 - 588*x2*x3^2*x4 - 350*x2*x3^2 + 8064*x2*x4^2 + 4800*x2*x4 + 5233*x3^4 - 10224*x3^3 + 289*x3^2*x4^2 - 1344*x3^2*x4 + 6768*x3^2*x5 + 5041*x3^2 - 850*x3*x4^2*x5 - 3094*x3*x4*x5^2 - 170*x3*x4*x5 - 6674*x3*x5 + 625*x4^2*x5^2 + 9216*x4^2 + 4550*x4*x5^3 + 250*x4*x5^2 + 8281*x5^4 + 910*x5^3 + 2234*x5^2
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
