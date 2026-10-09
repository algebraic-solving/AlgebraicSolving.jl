# Benchmark for `equidimensional_decomposition` on the "Sys2946" system.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AlgebraicSolving
using FileWatching.Pidfile: mkpidlock

const NAME = "Sys2946"

R, (a1, a2, b1, b2, c1, c2, o1, o2, x1, x2, y1, y2, z1, z2) = polynomial_ring(QQ, ["a1", "a2", "b1", "b2", "c1", "c2", "o1", "o2", "x1", "x2", "y1", "y2", "z1", "z2"])

F = [
    a1,
    a2,
    b1 - 1,
    b2,
    a1*b2 - a1*x2 - a2*b1 + a2*x1 + b1*x2 - b2*x1,
    b1*c2 - b1*y2 - b2*c1 + b2*y1 + c1*y2 - c2*y1,
    a1*c2 - a1*z2 - a2*c1 + a2*z1 + c1*z2 - c2*z1,
    c1*o2 - c1*x2 - c2*o1 + c2*x1 + o1*x2 - o2*x1,
    a1*o2 - a1*y2 - a2*o1 + a2*y1 + o1*y2 - o2*y1,
    b1*o2 - b1*z2 - b2*o1 + b2*z1 + o1*z2 - o2*z1
]

# Reducing the rational system is far cheaper than building the same
# polynomials again over the prime field.
prime = Int32(AlgebraicSolving.Nemo.rand_bits_prime(ZZ, 31))
Rp, _ = polynomial_ring(GF(prime), ["a1", "a2", "b1", "b2", "c1", "c2", "o1", "o2", "x1", "x2", "y1", "y2", "z1", "z2"])
Fp = [AlgebraicSolving.reduce_mod_p(f, Rp) for f in F]

# Decompose a tiny system first so that the timings below measure the
# computation rather than compilation. The solver specialises on the number of
# variables, so the warm up uses as many of them as the benchmark itself.
function warmup(K)
    _, w = polynomial_ring(K, ["w$(i)" for i in 1:14])
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
        println(io, join((NAME, 14, 10,
                          round(time_qq, digits = 6),
                          round(time_gf, digits = 6),
                          round(ratio, digits = 4)), ","))
    end
end
println("appended to $(csv)")
