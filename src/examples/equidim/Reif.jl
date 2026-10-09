# Benchmark for `equidimensional_decomposition` on the "Reif" system.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AlgebraicSolving
using FileWatching.Pidfile: mkpidlock

const NAME = "Reif"

R, (x16, x15, x14, x13, x12, x11, x10, x9, x8, x7, x6, x5, x4, x3, x2, x1) = polynomial_ring(QQ, ["x16", "x15", "x14", "x13", "x12", "x11", "x10", "x9", "x8", "x7", "x6", "x5", "x4", "x3", "x2", "x1"])

F = [
    -x14*x6 + x14*x5 - x13*x6 + x13*x4 + x6,
    -x16*x6 + x16*x5 - x15*x6 + x15*x4,
    -x14*x9 + x14*x8 - x13*x9 + x13*x7 + x9,
    -x16*x9 + x16*x8 - x15*x9 + x15*x7 - 1,
    -x14*x12 + x14*x11 - x13*x12 + x13*x10 + x12,
    -x16*x12 + x16*x11 - x15*x12 + x15*x10,
    -x14*x3 + x14*x2 - x13*x3 + x13*x1 + x3,
    -x16*x3 + x16*x2 - x15*x3 + x15*x1,
    -x14*x6*x3 + x14*x5*x2 - x13*x6*x3 + x13*x4*x1 + x6*x3 - 1,
    -x16*x6*x3 + x16*x5*x2 - x15*x6*x3 + x15*x4*x1,
    -x14*x9*x3 + x14*x8*x2 - x13*x9*x3 + x13*x7*x1 + x9*x3,
    -x16*x9*x3 + x16*x8*x2 - x15*x9*x3 + x15*x7*x1,
    -x14*x12*x3 + x14*x11*x2 - x13*x12*x3 + x13*x10*x1 + x12*x3,
    -x16*x12*x3 + x16*x11*x2 - x15*x12*x3 + x15*x10*x1 - 1
]

# Reducing the rational system is far cheaper than building the same
# polynomials again over the prime field.
prime = Int32(AlgebraicSolving.Nemo.rand_bits_prime(ZZ, 31))
Rp, _ = polynomial_ring(GF(prime), ["x16", "x15", "x14", "x13", "x12", "x11", "x10", "x9", "x8", "x7", "x6", "x5", "x4", "x3", "x2", "x1"])
Fp = [AlgebraicSolving.reduce_mod_p(f, Rp) for f in F]

# Decompose a tiny system first so that the timings below measure the
# computation rather than compilation. The solver specialises on the number of
# variables, so the warm up uses as many of them as the benchmark itself.
function warmup(K)
    _, w = polynomial_ring(K, ["w$(i)" for i in 1:16])
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
        println(io, join((NAME, 16, 14,
                          round(time_qq, digits = 6),
                          round(time_gf, digits = 6),
                          round(ratio, digits = 4)), ","))
    end
end
println("appended to $(csv)")
