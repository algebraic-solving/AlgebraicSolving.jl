# Benchmark for `equidimensional_decomposition` on the "pseudosing_deg_2_nvars_6" system.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AlgebraicSolving
using FileWatching.Pidfile: mkpidlock

const NAME = "pseudosing_deg_2_nvars_6"

R, (x1, x2, x3, y1, y2, y3) = polynomial_ring(QQ, ["x1", "x2", "x3", "y1", "y2", "y3"])

F = [
    x3 - y3,
    x2 - y2,
    -7*x1^2 + 22*x1*x2 - 55*x1*x3 - 94*x1 + 87*x2^2 - 56*x2*x3 - 62*x3^2 + 97*x3 - 73,
    -4*x1^2 - 83*x1*x2 - 10*x1*x3 + 62*x1 - 82*x2^2 + 80*x2*x3 - 44*x2 + 71*x3^2 - 17*x3 - 75,
    -7*y1^2 + 22*y1*y2 - 55*y1*y3 - 94*y1 + 87*y2^2 - 56*y2*y3 - 62*y3^2 + 97*y3 - 73,
    -4*y1^2 - 83*y1*y2 - 10*y1*y3 + 62*y1 - 82*y2^2 + 80*y2*y3 - 44*y2 + 71*y3^2 - 17*y3 - 75
]

# Reducing the rational system is far cheaper than building the same
# polynomials again over the prime field.
prime = Int32(AlgebraicSolving.Nemo.rand_bits_prime(ZZ, 31))
Rp, _ = polynomial_ring(GF(prime), ["x1", "x2", "x3", "y1", "y2", "y3"])
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
