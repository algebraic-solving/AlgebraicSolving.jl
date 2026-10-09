# Benchmark for `equidimensional_decomposition` on the "jorge2" system.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AlgebraicSolving
using FileWatching.Pidfile: mkpidlock

const NAME = "jorge2"

R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])

F = [
    73960*x^3 + 46428*x^2*y - 320426*x^2*z - 210018*x^2 - 88867*x*y^2 - 163934*x*y*z + 721747*x*y - 184389*x*z^2 + 416981*x*z - 111106*x - 62940*y^3 + 356381*y^2*z + 146898*y^2 + 32282*y*z^2 + 118097*y*z - 377082*y - 27183*z^3 + 116973*z^2 - 153504*z + 56580,
    3038*x^2*y - 3686*x^2*z + 2288*x*y^2 - 16544*x*y*z - 27166*x*y - 3344*x*z^2 + 4168*x*z + 315*y^3 + 4111*y^2*z - 2769*y^2 + 157*y*z^2 + 13942*y*z + 25806*y - 663*z^3 + 1527*z^2 - 690*z,
    650*x^2*y + 3350*x^2*z + 195*x*y^2 - 8450*x*y*z + 13390*x*y + 845*x*z^2 - 1630*x*z - 260*y^3 - 3410*y^2*z - 276*y^2 - 390*y*z^2 + 8372*y*z - 15088*y,
    -225*x^2*y + 1685*x^2*z - 705*x*y^2 - 450*x*y*z + 8056*x*y + 345*x*z^2 - 705*x*z - 420*y^3 + 1325*y^2*z + 918*y^2 + 645*y*z^2 + 1527*y*z - 5658*y
]

# Reducing the rational system is far cheaper than building the same
# polynomials again over the prime field.
prime = Int32(AlgebraicSolving.Nemo.rand_bits_prime(ZZ, 31))
Rp, _ = polynomial_ring(GF(prime), ["x", "y", "z"])
Fp = [AlgebraicSolving.reduce_mod_p(f, Rp) for f in F]

# Decompose a tiny system first so that the timings below measure the
# computation rather than compilation. The solver specialises on the number of
# variables, so the warm up uses as many of them as the benchmark itself.
function warmup(K)
    _, w = polynomial_ring(K, ["w$(i)" for i in 1:3])
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
        println(io, join((NAME, 3, 4,
                          round(time_qq, digits = 6),
                          round(time_gf, digits = 6),
                          round(ratio, digits = 4)), ","))
    end
end
println("appended to $(csv)")
