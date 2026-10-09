# Benchmark for `equidimensional_decomposition` on the "pseudosing_deg_2_nvars_10" system.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AlgebraicSolving
using FileWatching.Pidfile: mkpidlock

const NAME = "pseudosing_deg_2_nvars_10"

R, (x1, x2, x3, x4, x5, y1, y2, y3, y4, y5) = polynomial_ring(QQ, ["x1", "x2", "x3", "x4", "x5", "y1", "y2", "y3", "y4", "y5"])

F = [
    x5 - y5,
    x4 - y4,
    72*x1^2 - 87*x1*x2 + 47*x1*x3 - 90*x1*x4 + 43*x1*x5 + 92*x1 - 91*x2^2 - 88*x2*x3 - 48*x2*x4 + 53*x2*x5 - 28*x2 + 5*x3^2 + 13*x3*x4 - 10*x3*x5 - 82*x3 + 71*x4^2 + 16*x4*x5 + 83*x4 + 9*x5^2 - 60*x5 - 83,
    98*x1^2 - 48*x1*x2 - 19*x1*x3 + 62*x1*x4 + 37*x1*x5 + 5*x1 + 96*x2^2 - 17*x2*x3 + 25*x2*x4 + 91*x2*x5 + 98*x3^2 - 64*x3*x4 + 64*x3*x5 - 90*x3 - 60*x4^2 - 34*x4*x5 - 13*x4 + 44*x5^2 - 2*x5 + 71,
    -47*x1^2 - 39*x1*x2 - 53*x1*x3 - 72*x1*x4 - 97*x1*x5 + 33*x1 + 10*x2^2 + 7*x2*x3 - 89*x2*x4 + 65*x2*x5 + 12*x2 - 25*x3^2 - 96*x3*x4 + 50*x3*x5 - 60*x3 - 42*x4^2 + 7*x4*x5 - 89*x4 - 70*x5^2 + 34*x5 - 68,
    -60*x1^2 + 16*x1*x2 + 52*x1*x3 - 20*x1*x4 - 4*x1*x5 - 89*x1 - 77*x2^2 + 69*x2*x3 + 80*x2*x4 + 28*x2*x5 - 42*x2 - 33*x3^2 + 21*x3*x4 - 35*x3*x5 + 97*x3 + 30*x4^2 - 64*x4*x5 + 89*x4 - 16*x5^2 + 59*x5 - 69,
    72*y1^2 - 87*y1*y2 + 47*y1*y3 - 90*y1*y4 + 43*y1*y5 + 92*y1 - 91*y2^2 - 88*y2*y3 - 48*y2*y4 + 53*y2*y5 - 28*y2 + 5*y3^2 + 13*y3*y4 - 10*y3*y5 - 82*y3 + 71*y4^2 + 16*y4*y5 + 83*y4 + 9*y5^2 - 60*y5 - 83,
    98*y1^2 - 48*y1*y2 - 19*y1*y3 + 62*y1*y4 + 37*y1*y5 + 5*y1 + 96*y2^2 - 17*y2*y3 + 25*y2*y4 + 91*y2*y5 + 98*y3^2 - 64*y3*y4 + 64*y3*y5 - 90*y3 - 60*y4^2 - 34*y4*y5 - 13*y4 + 44*y5^2 - 2*y5 + 71,
    -47*y1^2 - 39*y1*y2 - 53*y1*y3 - 72*y1*y4 - 97*y1*y5 + 33*y1 + 10*y2^2 + 7*y2*y3 - 89*y2*y4 + 65*y2*y5 + 12*y2 - 25*y3^2 - 96*y3*y4 + 50*y3*y5 - 60*y3 - 42*y4^2 + 7*y4*y5 - 89*y4 - 70*y5^2 + 34*y5 - 68,
    -60*y1^2 + 16*y1*y2 + 52*y1*y3 - 20*y1*y4 - 4*y1*y5 - 89*y1 - 77*y2^2 + 69*y2*y3 + 80*y2*y4 + 28*y2*y5 - 42*y2 - 33*y3^2 + 21*y3*y4 - 35*y3*y5 + 97*y3 + 30*y4^2 - 64*y4*y5 + 89*y4 - 16*y5^2 + 59*y5 - 69
]

# Reducing the rational system is far cheaper than building the same
# polynomials again over the prime field.
prime = Int32(AlgebraicSolving.Nemo.rand_bits_prime(ZZ, 31))
Rp, _ = polynomial_ring(GF(prime), ["x1", "x2", "x3", "x4", "x5", "y1", "y2", "y3", "y4", "y5"])
Fp = [AlgebraicSolving.reduce_mod_p(f, Rp) for f in F]

# Decompose a tiny system first so that the timings below measure the
# computation rather than compilation. The solver specialises on the number of
# variables, so the warm up uses as many of them as the benchmark itself.
function warmup(K)
    _, w = polynomial_ring(K, ["w$(i)" for i in 1:10])
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
        println(io, join((NAME, 10, 10,
                          round(time_qq, digits = 6),
                          round(time_gf, digits = 6),
                          round(ratio, digits = 4)), ","))
    end
end
println("appended to $(csv)")
