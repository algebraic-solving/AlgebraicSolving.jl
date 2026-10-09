# Benchmark for `equidimensional_decomposition` on the "singhypcritpoints_deg_2_nvars_4" system.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AlgebraicSolving
using FileWatching.Pidfile: mkpidlock

const NAME = "singhypcritpoints_deg_2_nvars_4"

R, (x1, x2, x3, x4) = polynomial_ring(QQ, ["x1", "x2", "x3", "x4"])

F = [
    10816588*x1^3 + 38077162*x1^2*x2 + 25140479*x1^2*x3 - 124056256*x1^2*x4 - 9400295*x1^2 - 193635828*x1*x2^2 - 6790460*x1*x2*x3 + 250307870*x1*x2*x4 - 82452832*x1*x2 - 10319588*x1*x3^2 + 5839912*x1*x3*x4 + 22255617*x1*x3 + 112915639*x1*x4^2 + 115474646*x1*x4 - 17533608*x1 + 162291312*x2^3 - 126154506*x2^2*x3 - 52085058*x2^2*x4 + 157662837*x2^2 + 133237254*x2*x3^2 - 122040686*x2*x3*x4 - 2637098*x2*x3 - 132982128*x2*x4^2 - 183748338*x2*x4 + 40968512*x2 + 16298428*x3^3 + 46623254*x3^2*x4 + 119183502*x3^2 - 26720875*x3*x4^2 - 131012023*x3*x4 + 10206682*x3 - 16047619*x4^3 - 88321124*x4^2 - 58141098*x4 + 12230767,
    -19696652*x1^3 + 25140479*x1^2*x2 - 38705922*x1^2*x3 + 30086207*x1^2*x4 - 1289811*x1^2 - 3395230*x1*x2^2 - 20639176*x1*x2*x3 + 5839912*x1*x2*x4 + 22255617*x1*x2 + 36342096*x1*x3^2 + 64629228*x1*x3*x4 + 109167336*x1*x3 - 5646418*x1*x4^2 + 150465505*x1*x4 + 133696391*x1 - 42051502*x2^3 + 133237254*x2^2*x3 - 61020343*x2^2*x4 - 1318549*x2^2 + 48895284*x2*x3^2 + 93246508*x2*x3*x4 + 238367004*x2*x3 - 26720875*x2*x4^2 - 131012023*x2*x4 + 10206682*x2 + 60295600*x3^3 + 28934928*x3^2*x4 + 124816848*x3^2 - 19580172*x3*x4^2 + 116104228*x3*x4 + 145409424*x3 - 9372130*x4^3 - 70772721*x4^2 - 111938286*x4 - 23665371,
    12014836*x1^3 - 124056256*x1^2*x2 + 30086207*x1^2*x3 - 31723130*x1^2*x4 + 64756487*x1^2 + 125153935*x1*x2^2 + 5839912*x1*x2*x3 + 225831278*x1*x2*x4 + 115474646*x1*x2 + 32314614*x1*x3^2 - 11292836*x1*x3*x4 + 150465505*x1*x3 - 11471604*x1*x4^2 - 38448834*x1*x4 - 71771275*x1 - 17361686*x2^3 - 61020343*x2^2*x3 - 132982128*x2^2*x4 - 91874169*x2^2 + 46623254*x2*x3^2 - 53441750*x2*x3*x4 - 131012023*x2*x3 - 48142857*x2*x4^2 - 176642248*x2*x4 - 58141098*x2 + 9644976*x3^3 - 19580172*x3^2*x4 + 58052114*x3^2 - 28116390*x3*x4^2 - 141545442*x3*x4 - 111938286*x3 + 18646000*x4^3 + 9369657*x4^2 - 11893932*x4 + 7665457,
    2297232*x1^4 + 10816588*x1^3*x2 - 19696652*x1^3*x3 + 12014836*x1^3*x4 - 39813328*x1^3 + 19038581*x1^2*x2^2 + 25140479*x1^2*x2*x3 - 124056256*x1^2*x2*x4 - 9400295*x1^2*x2 - 19352961*x1^2*x3^2 + 30086207*x1^2*x3*x4 - 1289811*x1^2*x3 - 15861565*x1^2*x4^2 + 64756487*x1^2*x4 + 99968607*x1^2 - 64545276*x1*x2^3 - 3395230*x1*x2^2*x3 + 125153935*x1*x2^2*x4 - 41226416*x1*x2^2 - 10319588*x1*x2*x3^2 + 5839912*x1*x2*x3*x4 + 22255617*x1*x2*x3 + 112915639*x1*x2*x4^2 + 115474646*x1*x2*x4 - 17533608*x1*x2 + 12114032*x1*x3^3 + 32314614*x1*x3^2*x4 + 54583668*x1*x3^2 - 5646418*x1*x3*x4^2 + 150465505*x1*x3*x4 + 133696391*x1*x3 - 3823868*x1*x4^3 - 19224417*x1*x4^2 - 71771275*x1*x4 - 46338000*x1 + 40572828*x2^4 - 42051502*x2^3*x3 - 17361686*x2^3*x4 + 52554279*x2^3 + 66618627*x2^2*x3^2 - 61020343*x2^2*x3*x4 - 1318549*x2^2*x3 - 66491064*x2^2*x4^2 - 91874169*x2^2*x4 + 20484256*x2^2 + 16298428*x2*x3^3 + 46623254*x2*x3^2*x4 + 119183502*x2*x3^2 - 26720875*x2*x3*x4^2 - 131012023*x2*x3*x4 + 10206682*x2*x3 - 16047619*x2*x4^3 - 88321124*x2*x4^2 - 58141098*x2*x4 + 12230767*x2 + 15073900*x3^4 + 9644976*x3^3*x4 + 41605616*x3^3 - 9790086*x3^2*x4^2 + 58052114*x3^2*x4 + 72704712*x3^2 - 9372130*x3*x4^3 - 70772721*x3*x4^2 - 111938286*x3*x4 - 23665371*x3 + 4661500*x4^4 + 3123219*x4^3 - 5946966*x4^2 + 7665457*x4 + 9313394
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
