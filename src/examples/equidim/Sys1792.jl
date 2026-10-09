# Benchmark for `equidimensional_decomposition` on the "Sys1792" system.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AlgebraicSolving
using FileWatching.Pidfile: mkpidlock

const NAME = "Sys1792"

R, (x_1, x_2, x_3, x_4) = polynomial_ring(QQ, ["x_1", "x_2", "x_3", "x_4"])

F = [
    18*x_1^3*x_2*x_4 - 12*x_1^3*x_3^2 - 12*x_1^2*x_2*x_3 + 2*x_1^2*x_3*x_4^2 - 192*x_1*x_2^2 - 80*x_1*x_2*x_4^2 + 54*x_1*x_3^2*x_4 + 288*x_2*x_3*x_4 - 108*x_3^3 - 8*x_3*x_4^3,
    -108*x_1^3*x_2^2 + 54*x_1^2*x_2*x_3*x_4 - 12*x_1^2*x_3^3 + 288*x_1*x_2^2*x_4 - 12*x_1*x_2*x_3^2 - 8*x_1*x_2*x_4^3 + 2*x_1*x_3^2*x_4^2 - 192*x_2^2*x_3 - 80*x_2*x_3*x_4^2 + 18*x_3^3*x_4,
    18*x_1^3*x_2*x_3 + 144*x_1^2*x_2^2 - 12*x_1^2*x_2*x_4^2 + 2*x_1^2*x_3^2*x_4 - 160*x_1*x_2*x_3*x_4 + 18*x_1*x_3^3 - 256*x_2^2*x_4 + 144*x_2*x_3^2 + 64*x_2*x_4^3 - 12*x_3^2*x_4^2,
    -54*x_1^4*x_2 + 18*x_1^3*x_3*x_4 + 288*x_1^2*x_2*x_4 - 6*x_1^2*x_3^2 - 4*x_1^2*x_4^3 - 384*x_1*x_2*x_3 - 80*x_1*x_3*x_4^2 + 768*x_2^2 - 256*x_2*x_4^2 + 144*x_3^2*x_4 + 16*x_4^4,
    -27*x_1^4*x_2^2 + 18*x_1^3*x_2*x_3*x_4 - 4*x_1^3*x_3^3 + 144*x_1^2*x_2^2*x_4 - 6*x_1^2*x_2*x_3^2 - 4*x_1^2*x_2*x_4^3 + x_1^2*x_3^2*x_4^2 - 192*x_1*x_2^2*x_3 - 80*x_1*x_2*x_3*x_4^2 + 18*x_1*x_3^3*x_4 + 256*x_2^3 - 128*x_2^2*x_4^2 + 144*x_2*x_3^2*x_4 + 16*x_2*x_4^4 - 27*x_3^4 - 4*x_3^2*x_4^3
]

# Reducing the rational system is far cheaper than building the same
# polynomials again over the prime field.
prime = Int32(AlgebraicSolving.Nemo.rand_bits_prime(ZZ, 31))
Rp, _ = polynomial_ring(GF(prime), ["x_1", "x_2", "x_3", "x_4"])
Fp = [AlgebraicSolving.reduce_mod_p(f, Rp) for f in F]

# Decompose a tiny system first so that the timings below measure the
# computation rather than compilation. The solver specialises on the number of
# variables, so the warm up uses as many of them as the benchmark itself.
function warmup(K)
    _, w = polynomial_ring(K, ["w$(i)" for i in 1:4])
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
        println(io, join((NAME, 4, 5,
                          round(time_qq, digits = 6),
                          round(time_gf, digits = 6),
                          round(ratio, digits = 4)), ","))
    end
end
println("appended to $(csv)")
