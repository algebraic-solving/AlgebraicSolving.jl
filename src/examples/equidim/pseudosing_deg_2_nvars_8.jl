# Benchmark for `equidimensional_decomposition` on the "pseudosing_deg_2_nvars_8" system.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AlgebraicSolving
using FileWatching.Pidfile: mkpidlock

const NAME = "pseudosing_deg_2_nvars_8"

R, (x1, x2, x3, x4, y1, y2, y3, y4) = polynomial_ring(QQ, ["x1", "x2", "x3", "x4", "y1", "y2", "y3", "y4"])

F = [
    x4 - y4,
    x3 - y3,
    -10*x1^2 - 7*x1*x2 - 40*x1*x3 + 42*x1*x4 - 50*x1 + 23*x2^2 + 75*x2*x3 - 92*x2*x4 + 6*x2 + 74*x3^2 + 72*x3*x4 + 37*x3 - 23*x4^2 + 87*x4 + 44,
    29*x1^2 + 98*x1*x2 - 23*x1*x3 + 10*x1*x4 - 61*x1 - 8*x2^2 - 29*x2*x3 + 95*x2*x4 + 11*x2 - 49*x3^2 - 47*x3*x4 + 40*x3 - 81*x4^2 + 91*x4 + 68,
    -10*x1^2 + 31*x1*x2 - 51*x1*x3 + 77*x1*x4 + 95*x1 + x2^2 + x2*x3 + 55*x2*x4 - 28*x2 + 16*x3^2 + 30*x3*x4 - 27*x3 - 15*x4^2 - 59*x4 - 96,
    -10*y1^2 - 7*y1*y2 - 40*y1*y3 + 42*y1*y4 - 50*y1 + 23*y2^2 + 75*y2*y3 - 92*y2*y4 + 6*y2 + 74*y3^2 + 72*y3*y4 + 37*y3 - 23*y4^2 + 87*y4 + 44,
    29*y1^2 + 98*y1*y2 - 23*y1*y3 + 10*y1*y4 - 61*y1 - 8*y2^2 - 29*y2*y3 + 95*y2*y4 + 11*y2 - 49*y3^2 - 47*y3*y4 + 40*y3 - 81*y4^2 + 91*y4 + 68,
    -10*y1^2 + 31*y1*y2 - 51*y1*y3 + 77*y1*y4 + 95*y1 + y2^2 + y2*y3 + 55*y2*y4 - 28*y2 + 16*y3^2 + 30*y3*y4 - 27*y3 - 15*y4^2 - 59*y4 - 96
]

# Reducing the rational system is far cheaper than building the same
# polynomials again over the prime field.
prime = Int32(AlgebraicSolving.Nemo.rand_bits_prime(ZZ, 31))
Rp, _ = polynomial_ring(GF(prime), ["x1", "x2", "x3", "x4", "y1", "y2", "y3", "y4"])
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
