# Benchmark for `equidimensional_decomposition` on the "MontesS16" system.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AlgebraicSolving
using FileWatching.Pidfile: mkpidlock

const NAME = "MontesS16"

R, (w12, w13, w14, w15, w23, w25, w34, w45, w26, w36, w46, w56, x, y, z) = polynomial_ring(QQ, ["w12", "w13", "w14", "w15", "w23", "w25", "w34", "w45", "w26", "w36", "w46", "w56", "x", "y", "z"])

F = [
    w12 + w14,
    w12 + w13,
    w12 + w15,
    w12 + w23 + w25 - w26*x + w26,
    w12 + w25 - w26*y + w26,
    w12 + w23 - w26*z + w26,
    w23 + w34 + w36*x,
    w13 + w34 - w36*y + w36,
    w23 + w36*z,
    w14 + w34 + w45 - w46*x + w46,
    w34 + w46*y,
    w45 + w56*z,
    w15 + w45 - w56*z + w56,
    w26*x - w26 + w36*x + w46*x - w46 + w56*x,
    w26*y - w26 + w36*y - w36 + w46*y + w56*y,
    w26*z - w26 + w36*z + w46*z + w56*z - w56
]

# Reducing the rational system is far cheaper than building the same
# polynomials again over the prime field.
prime = Int32(AlgebraicSolving.Nemo.rand_bits_prime(ZZ, 31))
Rp, _ = polynomial_ring(GF(prime), ["w12", "w13", "w14", "w15", "w23", "w25", "w34", "w45", "w26", "w36", "w46", "w56", "x", "y", "z"])
Fp = [AlgebraicSolving.reduce_mod_p(f, Rp) for f in F]

# Decompose a tiny system first so that the timings below measure the
# computation rather than compilation. The solver specialises on the number of
# variables, so the warm up uses as many of them as the benchmark itself.
function warmup(K)
    _, w = polynomial_ring(K, ["w$(i)" for i in 1:15])
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
        println(io, join((NAME, 15, 16,
                          round(time_qq, digits = 6),
                          round(time_gf, digits = 6),
                          round(ratio, digits = 4)), ","))
    end
end
println("appended to $(csv)")
