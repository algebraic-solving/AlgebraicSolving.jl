# Benchmark for `equidimensional_decomposition` on the "Sys2945" system.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AlgebraicSolving
using FileWatching.Pidfile: mkpidlock

const NAME = "Sys2945"

R, (a, t, u, v, w, x, y, z) = polynomial_ring(QQ, ["a", "t", "u", "v", "w", "x", "y", "z"])

F = [
    -x^2 + y^2,
    a*x - u*w - u*x + v*y,
    -w^2 + x*y,
    -t*u + u*w + v*y - x*z,
    -a^2 + v^2,
    -a*w - a*y + u*w + u*y,
    -a*y + w*z,
    t*w - y^2,
    a*w - a*x + v*w - x*z,
    -a*y + t*u + u*y - x*z,
    t^2 - x*y,
    -a*v + a*z + u*v - u*z,
    -a*w + y*z,
    t*z - v*w - v*y + x*z,
    a*t - v*x,
    a^2 - z^2,
    -a*t + a*y - t*z + v*y,
    a*x - t*v,
    a*t + a*x - t*u - u*x,
    -t*w + x^2,
    -a^2 + v^2,
    -a*v + a*z + u*v - u*z,
    a*t - t*u - u*y + v*w,
    a^2 - z^2,
    -t*z - u*x + u*y + v*w,
    -a*w - t*z + u*w + u*x,
    -t^2 + w^2,
    a*u*y + u*v*x - w - x,
    a^2*x - a^2*y - a*u*x + a*u*y - a*v*x + a*y*z + u^2*x - u^2*y,
    a^2 + a*t*y + a*w*y + u*w*x - u*w*y - v^2 - v*x*y - v*y^2,
    -a*u*x + t - u*y*z + y,
    a*v - a*z + t*v*x - u*v + u*z - v*x*y - w*y*z + x*y*z,
    -a^2 - a*t*x - a*w*x + t*u*x - t*u*y + x^2*z + x*y*z + z^2,
    t*w*x - t*w*y + t*x*y + t*z - v*w - v*y - w*x*y + x*z,
    -u^2*v + x^2 + z,
    a*v*w + a*v*y + t - u*v*w - u*v*y - v^2*w - v^2*y + x,
    a^2*u - a*u^2 + a - u*v*z + w*y,
    a*u*y + a*v*w + a*v*x - u*v*w - v^2*x - v*w*z,
    -a*u*w + a*u*x - t*u*v + u*y*z,
    a^2 + a*w*x + a*w*y + t*u*y - t*v*w - u*w*y - v^2 - v*w^2,
    -t*w + u^2*v - z,
    a*t*v + a*u*w + a*v*y - t*v^2 - u*v*y - v*y*z,
    t*w - u^2*z + v,
    -a*t*u - a*w*z - a*x*z + u*x*z + v*x*z + w*z^2,
    a*v - a*z - t*v*w + t*v*x + t*w*z - u*v + u*z - w*y*z,
    u^2*v - x*y - z,
    -a^2*u + a*u^2 - a - t*x + u*v*z,
    a*u*w + t*u*v - t - y,
    -a^2*u + a*u^2 - a - t*x + u*v*z,
    u^2*z - v - y^2,
    a*t*u - a*u*y - u*v*x + u*w*z,
    -a*t*z - a*u*x - a*y*z + t*u*z + t*v*z + y*z^2,
    -a*t*z - a*x*z + t*u*z + t*z^2 + u*x*z - w + x*z^2 - y,
    -a^2 - a*t*x - a*t*y + t^2*z + t*u*x + t*w*z - u*w*x + z^2,
    t^2 - u^2*v + z,
    u^2*z - v - w^2,
    a^2*t - a^2*w - a*t*u - a*t*v + a*u*w + a*w*z + t*u^2 - u^2*w,
    a^2*u - a*u^2 + a - u*v*z + w*y,
    -u^2*z + v + x*y,
    -a*t*u - u*w*z + w + x
]

# Reducing the rational system is far cheaper than building the same
# polynomials again over the prime field.
prime = Int32(AlgebraicSolving.Nemo.rand_bits_prime(ZZ, 31))
Rp, _ = polynomial_ring(GF(prime), ["a", "t", "u", "v", "w", "x", "y", "z"])
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
        println(io, join((NAME, 8, 60,
                          round(time_qq, digits = 6),
                          round(time_gf, digits = 6),
                          round(ratio, digits = 4)), ","))
    end
end
println("appended to $(csv)")
