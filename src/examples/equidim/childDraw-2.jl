# Benchmark for `equidimensional_decomposition` on the "childDraw-2" system.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AlgebraicSolving
using FileWatching.Pidfile: mkpidlock

const NAME = "childDraw-2"

R, (a20, a32, a21, a31, a22, a30, a23, a35, a34, a33) = polynomial_ring(QQ, ["a20", "a32", "a21", "a31", "a22", "a30", "a23", "a35", "a34", "a33"])

F = [
    -80*a23 + 855*a35 + 180*a34,
    210*a35 - 210,
    16*a20*a32 + 18*a21*a31 + 20*a22*a30,
    7*a20*a31 + 8*a21*a30,
    40*a20*a34 + 48*a32*a22 + 44*a21*a33 + 52*a31*a23 + 280*a30,
    27*a20*a33 + 30*a32*a21 + 33*a31*a22 + 36*a30*a23,
    55*a20*a35 + 70*a32*a23 + 60*a21*a34 + 375*a31 + 65*a22*a33 + 80*a30,
    -170*a20 + 480*a32 + 78*a21*a35 + 102*a31 + 84*a22*a34 + 90*a23*a33,
    -114*a22 + 136*a23*a35 + 720*a34 + 152*a33,
    126*a32 - 144*a21 + 105*a22*a35 + 112*a23*a34 + 595*a33
]

# Reducing the rational system is far cheaper than building the same
# polynomials again over the prime field.
prime = Int32(AlgebraicSolving.Nemo.rand_bits_prime(ZZ, 31))
Rp, _ = polynomial_ring(GF(prime), ["a20", "a32", "a21", "a31", "a22", "a30", "a23", "a35", "a34", "a33"])
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
