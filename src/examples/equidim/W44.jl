# Benchmark for `equidimensional_decomposition` on the "W44" system.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AlgebraicSolving
using FileWatching.Pidfile: mkpidlock

const NAME = "W44"

R, (A0, A2, A3, A4, A5, B0, B1, B2, B3, B4, B5, C0, C1, C2, C3, C4, C5) = polynomial_ring(QQ, ["A0", "A2", "A3", "A4", "A5", "B0", "B1", "B2", "B3", "B4", "B5", "C0", "C1", "C2", "C3", "C4", "C5"])

F = [
    A2*B2,
    A4*B4,
    A5*B5,
    A2*B1 + B2,
    A4*B1 + B4,
    A2*B4 + A4*B2,
    A3*B5 + A5*B3,
    A4*B5 + A5*B4,
    A3*B5 + A5*B3 + 2*A5*B5,
    2*A3*B3 + A3*B5 + A5*B3,
    A3*B3 + A3*B5 + A5*B3 + A5*B5,
    A3*B4 + A4*B3 + A5*B1 + B5,
    A0*B2 + A2*B0 + A2*B1 + A2*B4 + A4*B2 + B2 + C2,
    A0*B1 + A2*B3 + A3*B2 + A4*B1 + B0 + 2*B1 + B4 + C1,
    A0*B4 + A2*B5 + A4*B0 + A4*B1 + 2*A4*B4 + A5*B2 + B4 + C4,
    A0*B3 + A3*B0 + 2*A3*B1 + A3*B4 + A4*B3 + A5*B1 + 2*B3 + B5 + C3,
    A0*B5 + A3*B4 + A4*B3 + 2*A4*B5 + A5*B0 + A5*B1 + 2*A5*B4 + B5 + C5,
    A0*B3 + A0*B5 + A3*B0 + A3*B1 + A3*B4 + A4*B3 + A4*B5 + A5*B0 + A5*B1 + A5*B4 + B3 + B5 + C3 + C5 - 1,
    A0*B0 + A0*B1 + A0*B4 + A2*B3 + A2*B5 + A3*B2 + A4*B0 + A4*B1 + A4*B4 + A5*B2 + B0 + B1 + B4 + C0 + C1 + C4
]

# Reducing the rational system is far cheaper than building the same
# polynomials again over the prime field.
prime = Int32(AlgebraicSolving.Nemo.rand_bits_prime(ZZ, 31))
Rp, _ = polynomial_ring(GF(prime), ["A0", "A2", "A3", "A4", "A5", "B0", "B1", "B2", "B3", "B4", "B5", "C0", "C1", "C2", "C3", "C4", "C5"])
Fp = [AlgebraicSolving.reduce_mod_p(f, Rp) for f in F]

# Decompose a tiny system first so that the timings below measure the
# computation rather than compilation. The solver specialises on the number of
# variables, so the warm up uses as many of them as the benchmark itself.
function warmup(K)
    _, w = polynomial_ring(K, ["w$(i)" for i in 1:17])
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
        println(io, join((NAME, 17, 19,
                          round(time_qq, digits = 6),
                          round(time_gf, digits = 6),
                          round(ratio, digits = 4)), ","))
    end
end
println("appended to $(csv)")
