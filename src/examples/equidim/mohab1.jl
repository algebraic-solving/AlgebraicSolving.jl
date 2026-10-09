# Benchmark for `equidimensional_decomposition` on the "mohab1" system.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AlgebraicSolving
using FileWatching.Pidfile: mkpidlock

const NAME = "mohab1"

R, (u, v, w, x) = polynomial_ring(QQ, ["u", "v", "w", "x"])

F = [
    160*u^4*v^3*w^3 - 64*u^3*v^4*w^4*x - 96*u^3*v^4*w^2*x - 192*u^3*v^4*w^2 + 108*u^3*v^4*x - 96*u^3*v^2*w^4*x - 192*u^3*v^2*w^4 - 24*u^3*v^2*w^2*x + 108*u^3*w^4*x + 48*u^2*v^5*w^3*x^2 + 48*u^2*v^5*w^3*x - 36*u^2*v^5*w*x^2 - 36*u^2*v^5*w*x + 72*u^2*v^5*w + 48*u^2*v^3*w^5*x^2 + 48*u^2*v^3*w^5*x + 48*u^2*v^3*w^3*x^2 + 456*u^2*v^3*w^3*x + 48*u^2*v^3*w^3 - 72*u^2*v^3*w*x^2 - 144*u^2*v^3*w*x - 36*u^2*v*w^5*x^2 - 36*u^2*v*w^5*x + 72*u^2*v*w^5 - 72*u^2*v*w^3*x^2 - 144*u^2*v*w^3*x + 72*u^2*v*w*x^2 - 8*u*v^6*w^2*x^3 + 16*u*v^6*w^2*x^2 - 8*u*v^6*w^2*x + 8*u*v^6*x^3 - 24*u*v^6*x^2 + 24*u*v^6*x - 8*u*v^6 - 32*u*v^4*w^4*x^3 - 128*u*v^4*w^4*x^2 + 16*u*v^4*w^4*x + 36*u*v^4*w^2*x^3 - 120*u*v^4*w^2*x + 72*u*v^4*w^2 - 12*u*v^4*x^3 + 96*u*v^4*x^2 - 84*u*v^4*x - 8*u*v^2*w^6*x^3 + 16*u*v^2*w^6*x^2 - 8*u*v^2*w^6*x + 36*u*v^2*w^4*x^3 - 120*u*v^2*w^4*x + 72*u*v^2*w^4 + 168*u*v^2*w^2*x - 12*u*v^2*x^3 - 24*u*v^2*x^2 + 8*u*w^6*x^3 - 24*u*w^6*x^2 + 24*u*w^6*x - 8*u*w^6 - 12*u*w^4*x^3 + 96*u*w^4*x^2 - 84*u*w^4*x - 12*u*w^2*x^3 - 24*u*w^2*x^2 + 8*u*x^3 + 4*v^5*w^3*x^4 + 12*v^5*w^3*x^3 - 16*v^5*w^3*x - 4*v^5*w*x^4 - 4*v^5*w*x^3 - 4*v^5*w*x^2 + 36*v^5*w*x - 24*v^5*w + 4*v^3*w^5*x^4 + 12*v^3*w^5*x^3 - 16*v^3*w^5*x - 16*v^3*w^3*x^4 - 16*v^3*w^3*x^3 + 80*v^3*w^3*x^2 - 72*v^3*w^3*x - 48*v^3*w^3 + 12*v^3*w*x^4 - 12*v^3*w*x^3 - 24*v^3*w*x^2 + 48*v^3*w*x - 4*v*w^5*x^4 - 4*v*w^5*x^3 - 4*v*w^5*x^2 + 36*v*w^5*x - 24*v*w^5 + 12*v*w^3*x^4 - 12*v*w^3*x^3 - 24*v*w^3*x^2 + 48*v*w^3*x - 8*v*w*x^4 + 24*v*w*x^3 - 24*v*w*x^2,
    96*u^5*v^2*w^3 - 64*u^4*v^3*w^4*x - 96*u^4*v^3*w^2*x - 192*u^4*v^3*w^2 + 108*u^4*v^3*x - 48*u^4*v*w^4*x - 96*u^4*v*w^4 - 12*u^4*v*w^2*x + 80*u^3*v^4*w^3*x^2 + 80*u^3*v^4*w^3*x - 60*u^3*v^4*w*x^2 - 60*u^3*v^4*w*x + 120*u^3*v^4*w + 48*u^3*v^2*w^5*x^2 + 48*u^3*v^2*w^5*x + 48*u^3*v^2*w^3*x^2 + 456*u^3*v^2*w^3*x + 48*u^3*v^2*w^3 - 72*u^3*v^2*w*x^2 - 144*u^3*v^2*w*x - 12*u^3*w^5*x^2 - 12*u^3*w^5*x + 24*u^3*w^5 - 24*u^3*w^3*x^2 - 48*u^3*w^3*x + 24*u^3*w*x^2 - 24*u^2*v^5*w^2*x^3 + 48*u^2*v^5*w^2*x^2 - 24*u^2*v^5*w^2*x + 24*u^2*v^5*x^3 - 72*u^2*v^5*x^2 + 72*u^2*v^5*x - 24*u^2*v^5 - 64*u^2*v^3*w^4*x^3 - 256*u^2*v^3*w^4*x^2 + 32*u^2*v^3*w^4*x + 72*u^2*v^3*w^2*x^3 - 240*u^2*v^3*w^2*x + 144*u^2*v^3*w^2 - 24*u^2*v^3*x^3 + 192*u^2*v^3*x^2 - 168*u^2*v^3*x - 8*u^2*v*w^6*x^3 + 16*u^2*v*w^6*x^2 - 8*u^2*v*w^6*x + 36*u^2*v*w^4*x^3 - 120*u^2*v*w^4*x + 72*u^2*v*w^4 + 168*u^2*v*w^2*x - 12*u^2*v*x^3 - 24*u^2*v*x^2 + 20*u*v^4*w^3*x^4 + 60*u*v^4*w^3*x^3 - 80*u*v^4*w^3*x - 20*u*v^4*w*x^4 - 20*u*v^4*w*x^3 - 20*u*v^4*w*x^2 + 180*u*v^4*w*x - 120*u*v^4*w + 12*u*v^2*w^5*x^4 + 36*u*v^2*w^5*x^3 - 48*u*v^2*w^5*x - 48*u*v^2*w^3*x^4 - 48*u*v^2*w^3*x^3 + 240*u*v^2*w^3*x^2 - 216*u*v^2*w^3*x - 144*u*v^2*w^3 + 36*u*v^2*w*x^4 - 36*u*v^2*w*x^3 - 72*u*v^2*w*x^2 + 144*u*v^2*w*x - 4*u*w^5*x^4 - 4*u*w^5*x^3 - 4*u*w^5*x^2 + 36*u*w^5*x - 24*u*w^5 + 12*u*w^3*x^4 - 12*u*w^3*x^3 - 24*u*w^3*x^2 + 48*u*w^3*x - 8*u*w*x^4 + 24*u*w*x^3 - 24*u*w*x^2 - 24*v^5*w^2*x^4 + 72*v^5*w^2*x^3 - 72*v^5*w^2*x^2 + 24*v^5*w^2*x + 24*v^5*x^4 - 96*v^5*x^3 + 144*v^5*x^2 - 96*v^5*x + 24*v^5 - 4*v^3*w^4*x^5 + 16*v^3*w^4*x^4 - 96*v^3*w^4*x^3 + 160*v^3*w^4*x^2 + 32*v^3*w^4*x + 8*v^3*w^2*x^5 + 24*v^3*w^2*x^3 - 176*v^3*w^2*x^2 + 96*v^3*w^2*x + 48*v^3*w^2 - 4*v^3*x^5 - 16*v^3*x^4 + 104*v^3*x^3 - 144*v^3*x^2 + 60*v^3*x - 8*v*w^6*x^4 + 24*v*w^6*x^3 - 24*v*w^6*x^2 + 8*v*w^6*x + 4*v*w^4*x^5 + 12*v*w^4*x^3 - 88*v*w^4*x^2 + 48*v*w^4*x + 24*v*w^4 - 8*v*w^2*x^5 + 16*v*w^2*x^4 - 32*v*w^2*x^3 + 144*v*w^2*x^2 - 156*v*w^2*x + 4*v*x^5 - 8*v*x^4 - 12*v*x^3 + 24*v*x^2,
    96*u^5*v^3*w^2 - 64*u^4*v^4*w^3*x - 48*u^4*v^4*w*x - 96*u^4*v^4*w - 96*u^4*v^2*w^3*x - 192*u^4*v^2*w^3 - 12*u^4*v^2*w*x + 108*u^4*w^3*x + 48*u^3*v^5*w^2*x^2 + 48*u^3*v^5*w^2*x - 12*u^3*v^5*x^2 - 12*u^3*v^5*x + 24*u^3*v^5 + 80*u^3*v^3*w^4*x^2 + 80*u^3*v^3*w^4*x + 48*u^3*v^3*w^2*x^2 + 456*u^3*v^3*w^2*x + 48*u^3*v^3*w^2 - 24*u^3*v^3*x^2 - 48*u^3*v^3*x - 60*u^3*v*w^4*x^2 - 60*u^3*v*w^4*x + 120*u^3*v*w^4 - 72*u^3*v*w^2*x^2 - 144*u^3*v*w^2*x + 24*u^3*v*x^2 - 8*u^2*v^6*w*x^3 + 16*u^2*v^6*w*x^2 - 8*u^2*v^6*w*x - 64*u^2*v^4*w^3*x^3 - 256*u^2*v^4*w^3*x^2 + 32*u^2*v^4*w^3*x + 36*u^2*v^4*w*x^3 - 120*u^2*v^4*w*x + 72*u^2*v^4*w - 24*u^2*v^2*w^5*x^3 + 48*u^2*v^2*w^5*x^2 - 24*u^2*v^2*w^5*x + 72*u^2*v^2*w^3*x^3 - 240*u^2*v^2*w^3*x + 144*u^2*v^2*w^3 + 168*u^2*v^2*w*x + 24*u^2*w^5*x^3 - 72*u^2*w^5*x^2 + 72*u^2*w^5*x - 24*u^2*w^5 - 24*u^2*w^3*x^3 + 192*u^2*w^3*x^2 - 168*u^2*w^3*x - 12*u^2*w*x^3 - 24*u^2*w*x^2 + 12*u*v^5*w^2*x^4 + 36*u*v^5*w^2*x^3 - 48*u*v^5*w^2*x - 4*u*v^5*x^4 - 4*u*v^5*x^3 - 4*u*v^5*x^2 + 36*u*v^5*x - 24*u*v^5 + 20*u*v^3*w^4*x^4 + 60*u*v^3*w^4*x^3 - 80*u*v^3*w^4*x - 48*u*v^3*w^2*x^4 - 48*u*v^3*w^2*x^3 + 240*u*v^3*w^2*x^2 - 216*u*v^3*w^2*x - 144*u*v^3*w^2 + 12*u*v^3*x^4 - 12*u*v^3*x^3 - 24*u*v^3*x^2 + 48*u*v^3*x - 20*u*v*w^4*x^4 - 20*u*v*w^4*x^3 - 20*u*v*w^4*x^2 + 180*u*v*w^4*x - 120*u*v*w^4 + 36*u*v*w^2*x^4 - 36*u*v*w^2*x^3 - 72*u*v*w^2*x^2 + 144*u*v*w^2*x - 8*u*v*x^4 + 24*u*v*x^3 - 24*u*v*x^2 - 8*v^6*w*x^4 + 24*v^6*w*x^3 - 24*v^6*w*x^2 + 8*v^6*w*x - 4*v^4*w^3*x^5 + 16*v^4*w^3*x^4 - 96*v^4*w^3*x^3 + 160*v^4*w^3*x^2 + 32*v^4*w^3*x + 4*v^4*w*x^5 + 12*v^4*w*x^3 - 88*v^4*w*x^2 + 48*v^4*w*x + 24*v^4*w - 24*v^2*w^5*x^4 + 72*v^2*w^5*x^3 - 72*v^2*w^5*x^2 + 24*v^2*w^5*x + 8*v^2*w^3*x^5 + 24*v^2*w^3*x^3 - 176*v^2*w^3*x^2 + 96*v^2*w^3*x + 48*v^2*w^3 - 8*v^2*w*x^5 + 16*v^2*w*x^4 - 32*v^2*w*x^3 + 144*v^2*w*x^2 - 156*v^2*w*x + 24*w^5*x^4 - 96*w^5*x^3 + 144*w^5*x^2 - 96*w^5*x + 24*w^5 - 4*w^3*x^5 - 16*w^3*x^4 + 104*w^3*x^3 - 144*w^3*x^2 + 60*w^3*x + 4*w*x^5 - 8*w*x^4 - 12*w*x^3 + 24*w*x^2,
    -16*u^4*v^4*w^4 - 24*u^4*v^4*w^2 + 27*u^4*v^4 - 24*u^4*v^2*w^4 - 6*u^4*v^2*w^2 + 27*u^4*w^4 + 32*u^3*v^5*w^3*x + 16*u^3*v^5*w^3 - 24*u^3*v^5*w*x - 12*u^3*v^5*w + 32*u^3*v^3*w^5*x + 16*u^3*v^3*w^5 + 32*u^3*v^3*w^3*x + 152*u^3*v^3*w^3 - 48*u^3*v^3*w*x - 48*u^3*v^3*w - 24*u^3*v*w^5*x - 12*u^3*v*w^5 - 48*u^3*v*w^3*x - 48*u^3*v*w^3 + 48*u^3*v*w*x - 12*u^2*v^6*w^2*x^2 + 16*u^2*v^6*w^2*x - 4*u^2*v^6*w^2 + 12*u^2*v^6*x^2 - 24*u^2*v^6*x + 12*u^2*v^6 - 48*u^2*v^4*w^4*x^2 - 128*u^2*v^4*w^4*x + 8*u^2*v^4*w^4 + 54*u^2*v^4*w^2*x^2 - 60*u^2*v^4*w^2 - 18*u^2*v^4*x^2 + 96*u^2*v^4*x - 42*u^2*v^4 - 12*u^2*v^2*w^6*x^2 + 16*u^2*v^2*w^6*x - 4*u^2*v^2*w^6 + 54*u^2*v^2*w^4*x^2 - 60*u^2*v^2*w^4 + 84*u^2*v^2*w^2 - 18*u^2*v^2*x^2 - 24*u^2*v^2*x + 12*u^2*w^6*x^2 - 24*u^2*w^6*x + 12*u^2*w^6 - 18*u^2*w^4*x^2 + 96*u^2*w^4*x - 42*u^2*w^4 - 18*u^2*w^2*x^2 - 24*u^2*w^2*x + 12*u^2*x^2 + 16*u*v^5*w^3*x^3 + 36*u*v^5*w^3*x^2 - 16*u*v^5*w^3 - 16*u*v^5*w*x^3 - 12*u*v^5*w*x^2 - 8*u*v^5*w*x + 36*u*v^5*w + 16*u*v^3*w^5*x^3 + 36*u*v^3*w^5*x^2 - 16*u*v^3*w^5 - 64*u*v^3*w^3*x^3 - 48*u*v^3*w^3*x^2 + 160*u*v^3*w^3*x - 72*u*v^3*w^3 + 48*u*v^3*w*x^3 - 36*u*v^3*w*x^2 - 48*u*v^3*w*x + 48*u*v^3*w - 16*u*v*w^5*x^3 - 12*u*v*w^5*x^2 - 8*u*v*w^5*x + 36*u*v*w^5 + 48*u*v*w^3*x^3 - 36*u*v*w^3*x^2 - 48*u*v*w^3*x + 48*u*v*w^3 - 32*u*v*w*x^3 + 72*u*v*w*x^2 - 48*u*v*w*x - 16*v^6*w^2*x^3 + 36*v^6*w^2*x^2 - 24*v^6*w^2*x + 4*v^6*w^2 + 16*v^6*x^3 - 48*v^6*x^2 + 48*v^6*x - 16*v^6 - 5*v^4*w^4*x^4 + 16*v^4*w^4*x^3 - 72*v^4*w^4*x^2 + 80*v^4*w^4*x + 8*v^4*w^4 + 10*v^4*w^2*x^4 + 18*v^4*w^2*x^2 - 88*v^4*w^2*x + 24*v^4*w^2 - 5*v^4*x^4 - 16*v^4*x^3 + 78*v^4*x^2 - 72*v^4*x + 15*v^4 - 16*v^2*w^6*x^3 + 36*v^2*w^6*x^2 - 24*v^2*w^6*x + 4*v^2*w^6 + 10*v^2*w^4*x^4 + 18*v^2*w^4*x^2 - 88*v^2*w^4*x + 24*v^2*w^4 - 20*v^2*w^2*x^4 + 32*v^2*w^2*x^3 - 48*v^2*w^2*x^2 + 144*v^2*w^2*x - 78*v^2*w^2 + 10*v^2*x^4 - 16*v^2*x^3 - 18*v^2*x^2 + 24*v^2*x + 16*w^6*x^3 - 48*w^6*x^2 + 48*w^6*x - 16*w^6 - 5*w^4*x^4 - 16*w^4*x^3 + 78*w^4*x^2 - 72*w^4*x + 15*w^4 + 10*w^2*x^4 - 16*w^2*x^3 - 18*w^2*x^2 + 24*w^2*x - 5*x^4 + 16*x^3 - 12*x^2
]

# Reducing the rational system is far cheaper than building the same
# polynomials again over the prime field.
prime = Int32(AlgebraicSolving.Nemo.rand_bits_prime(ZZ, 31))
Rp, _ = polynomial_ring(GF(prime), ["u", "v", "w", "x"])
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
        println(io, join((NAME, 4, 4,
                          round(time_qq, digits = 6),
                          round(time_gf, digits = 6),
                          round(ratio, digits = 4)), ","))
    end
end
println("appended to $(csv)")
