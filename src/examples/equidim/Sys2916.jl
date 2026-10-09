# Benchmark for `equidimensional_decomposition` on the "Sys2916" system.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AlgebraicSolving
using FileWatching.Pidfile: mkpidlock

const NAME = "Sys2916"

R, (A21, A31, A32, A41, A42, A43, B1, B2, B3, B4, C2, C3, C4) = polynomial_ring(QQ, ["A21", "A31", "A32", "A41", "A42", "A43", "B1", "B2", "B3", "B4", "C2", "C3", "C4"])

F = [
    B1 + B2 + B3 + B4 - 1,
    -A21 + C2,
    -A31 - A32 + C3,
    -A41 - A42 - A43 + C4,
    2*B2*C2 + 2*B3*C3 + 2*B4*C4 - 1,
    3*B2*C2^2 + 3*B3*C3^2 + 3*B4*C4^2 - 1,
    6*A32*B3*C2 + 6*A42*B4*C2 + 6*A43*B4*C3 - 1,
    4*B2*C2^3 + 4*B3*C3^3 + 4*B4*C4^3 - 1,
    8*A32*B3*C2*C3 + 8*A42*B4*C2*C4 + 8*A43*B4*C3*C4 - 1,
    12*A32*B3*C2^2 + 12*A42*B4*C2^2 + 12*A43*B4*C3^2 - 1,
    24*A32*A43*B4*C2 - 1
]

# Reducing the rational system is far cheaper than building the same
# polynomials again over the prime field.
prime = Int32(AlgebraicSolving.Nemo.rand_bits_prime(ZZ, 31))
Rp, _ = polynomial_ring(GF(prime), ["A21", "A31", "A32", "A41", "A42", "A43", "B1", "B2", "B3", "B4", "C2", "C3", "C4"])
Fp = [AlgebraicSolving.reduce_mod_p(f, Rp) for f in F]

# Decompose a tiny system first so that the timings below measure the
# computation rather than compilation. The solver specialises on the number of
# variables, so the warm up uses as many of them as the benchmark itself.
function warmup(K)
    _, w = polynomial_ring(K, ["w$(i)" for i in 1:13])
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
        println(io, join((NAME, 13, 11,
                          round(time_qq, digits = 6),
                          round(time_gf, digits = 6),
                          round(ratio, digits = 4)), ","))
    end
end
println("appended to $(csv)")
