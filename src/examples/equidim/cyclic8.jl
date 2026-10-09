# Benchmark for `equidimensional_decomposition` on the "cyclic8" system.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AlgebraicSolving
using FileWatching.Pidfile: mkpidlock

const NAME = "cyclic8"

R, (x1, x2, x3, x4, x5, x6, x7, x8) = polynomial_ring(QQ, ["x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8"])

F = [
    x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8,
    x1*x2 + x1*x8 + x2*x3 + x3*x4 + x4*x5 + x5*x6 + x6*x7 + x7*x8,
    x1*x2*x3 + x1*x2*x8 + x1*x7*x8 + x2*x3*x4 + x3*x4*x5 + x4*x5*x6 + x5*x6*x7 + x6*x7*x8,
    x1*x2*x3*x4 + x1*x2*x3*x8 + x1*x2*x7*x8 + x1*x6*x7*x8 + x2*x3*x4*x5 + x3*x4*x5*x6 + x4*x5*x6*x7 + x5*x6*x7*x8,
    x1*x2*x3*x4*x5 + x1*x2*x3*x4*x8 + x1*x2*x3*x7*x8 + x1*x2*x6*x7*x8 + x1*x5*x6*x7*x8 + x2*x3*x4*x5*x6 + x3*x4*x5*x6*x7 + x4*x5*x6*x7*x8,
    x1*x2*x3*x4*x5*x6 + x1*x2*x3*x4*x5*x8 + x1*x2*x3*x4*x7*x8 + x1*x2*x3*x6*x7*x8 + x1*x2*x5*x6*x7*x8 + x1*x4*x5*x6*x7*x8 + x2*x3*x4*x5*x6*x7 + x3*x4*x5*x6*x7*x8,
    x1*x2*x3*x4*x5*x6*x7 + x1*x2*x3*x4*x5*x6*x8 + x1*x2*x3*x4*x5*x7*x8 + x1*x2*x3*x4*x6*x7*x8 + x1*x2*x3*x5*x6*x7*x8 + x1*x2*x4*x5*x6*x7*x8 + x1*x3*x4*x5*x6*x7*x8 + x2*x3*x4*x5*x6*x7*x8,
    x1*x2*x3*x4*x5*x6*x7*x8 - 1
]

# Reducing the rational system is far cheaper than building the same
# polynomials again over the prime field.
prime = Int32(AlgebraicSolving.Nemo.rand_bits_prime(ZZ, 31))
Rp, _ = polynomial_ring(GF(prime), ["x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8"])
Fp = [AlgebraicSolving.reduce_mod_p(f, Rp) for f in F]

# Decompose a tiny system first so that the timings below measure the
# computation rather than compilation. The solver specialises on the number of
# variables, so the warm up uses as many of them as the benchmark itself.
function warmup(K)
    _, w = polynomial_ring(K, ["w$(i)" for i in 1:8])
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
        println(io, join((NAME, 8, 8,
                          round(time_qq, digits = 6),
                          round(time_gf, digits = 6),
                          round(ratio, digits = 4)), ","))
    end
end
println("appended to $(csv)")
