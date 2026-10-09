# Benchmark for `equidimensional_decomposition` on the "dgp6" system.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AlgebraicSolving
using FileWatching.Pidfile: mkpidlock

const NAME = "dgp6"

R, (s, q, p, o, n, m, l, k, j, h, g, f, e, d, c, b, a) = polynomial_ring(QQ, ["s", "q", "p", "o", "n", "m", "l", "k", "j", "h", "g", "f", "e", "d", "c", "b", "a"])

F = [
    g*a,
    q + p*n + m*a + j*g,
    l*b,
    q*n,
    p*l + o*l + l*a + k*b + g*b + c + b,
    q*l + p*k + p*g + p + o*k + o*g + o + n*b + m*b + l*j + k*a + h + g*a + f + d + a,
    s + q*k + q*g + q + p*n + p*m + o*n + o*m + n*a + m*a + k*j + j*g + j + e - 1,
    q*n + q*m + n*j + m*j,
    2*q*n + q*m + n*j,
    s + q*k + 2*q*g + q + p*n + o*n + 2*n*a + m*a + j*g,
    q*l + p*g + o*g + n*b + k*a + 2*g*a + d + a,
    l*a + g*b,
    q*g + n*a,
    q*m + n*j + 2*m*j,
    q + p*n + 2*p*m + o*m + m*a + k*j + j*g + 2*j + e,
    p*k + p*g + 2*p + o + m*b + l*j + f + a,
    p*l + b,
    q*m + n*j,
    p*g + a
]

# Reducing the rational system is far cheaper than building the same
# polynomials again over the prime field.
prime = Int32(AlgebraicSolving.Nemo.rand_bits_prime(ZZ, 31))
Rp, _ = polynomial_ring(GF(prime), ["s", "q", "p", "o", "n", "m", "l", "k", "j", "h", "g", "f", "e", "d", "c", "b", "a"])
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
