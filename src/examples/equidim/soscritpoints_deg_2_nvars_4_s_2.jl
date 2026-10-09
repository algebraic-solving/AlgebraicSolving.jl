# Benchmark for `equidimensional_decomposition` on the "soscritpoints_deg_2_nvars_4_s_2" system.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AlgebraicSolving
using FileWatching.Pidfile: mkpidlock

const NAME = "soscritpoints_deg_2_nvars_4_s_2"

R, (x1, x2, x3, x4) = polynomial_ring(QQ, ["x1", "x2", "x3", "x4"])

F = [
    13888*x1^2*x2 - 4828*x1*x2*x3 - 27016*x1*x2 - 10650*x1*x3^2 - 12000*x1*x3 + 16532*x2^3 + 7650*x2^2*x3 - 12806*x2*x3^2 + 18104*x2*x4^2 + 2992*x2*x4 + 5576*x2 + 6600*x3*x4 + 12300*x3,
    -11646*x1^2*x3 + 11360*x1^2 - 2414*x1*x2^2 - 21300*x1*x2*x3 - 12000*x1*x2 + 33756*x1*x3 - 6248*x1*x4 - 11644*x1 + 2550*x2^3 - 12806*x2^2*x3 + 6600*x2*x4 + 12300*x2 + 37636*x3^3 - 28324*x3*x4^2,
    16352*x1^2*x4 - 6248*x1*x3 - 25404*x1*x4 - 7040*x1 + 18104*x2^2*x4 + 1496*x2^2 + 6600*x2*x3 - 28324*x3^2*x4 + 21316*x4^3 + 3872*x4 + 7216,
    3136*x1^4 - 9744*x1^3 + 6944*x1^2*x2^2 - 5823*x1^2*x3^2 + 11360*x1^2*x3 + 8176*x1^2*x4^2 + 13969*x1^2 - 2414*x1*x2^2*x3 - 13508*x1*x2^2 - 10650*x1*x2*x3^2 - 12000*x1*x2*x3 + 16878*x1*x3^2 - 6248*x1*x3*x4 - 11644*x1*x3 - 12702*x1*x4^2 - 7040*x1*x4 - 13120*x1 + 4133*x2^4 + 2550*x2^3*x3 - 6403*x2^2*x3^2 + 9052*x2^2*x4^2 + 1496*x2^2*x4 + 2788*x2^2 + 6600*x2*x3*x4 + 12300*x2*x3 + 9409*x3^4 - 14162*x3^2*x4^2 + 5329*x4^4 + 1936*x4^2 + 7216*x4 + 6724
]

# Reducing the rational system is far cheaper than building the same
# polynomials again over the prime field.
prime = Int32(AlgebraicSolving.Nemo.rand_bits_prime(ZZ, 31))
Rp, _ = polynomial_ring(GF(prime), ["x1", "x2", "x3", "x4"])
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
