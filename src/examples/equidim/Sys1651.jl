# Benchmark for `equidimensional_decomposition` on the "Sys1651" system.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AlgebraicSolving
using FileWatching.Pidfile: mkpidlock

const NAME = "Sys1651"

R, (x_1, x_2, x_3) = polynomial_ring(QQ, ["x_1", "x_2", "x_3"])

F = [
    -ZZ("11208668022675968779547170612512235210349047243552000000000000")*x_1^7 - ZZ("31449700298556793984207504684027195888597735889737780000000000")*x_1^6 + ZZ("21816921759249384025948630953428245512680770541232190225000000")*x_1^5 + ZZ("45569733208624126847447964943715863629520914115150822264750000")*x_1^4 + ZZ("88805846299976113734136773560701804860459385122615385374907500")*x_1^3 - ZZ("86214424776973426748434098131092218366969015035499285477132575")*x_1^2 - ZZ("14718812231543448446504246005049744981109245109830198914473480")*x_1 + ZZ("42090682858706655449898213849821131604406527035405370572360500")*x_2 - ZZ("6484094511127944287766737675606481544721165689917908276697408"),
    -ZZ("460925738134370922105067572622357130113336819044908000000000000")*x_1^7 - ZZ("2008196121863672147132162589054374608823802723131731932500000000")*x_1^6 - ZZ("2299659624342332747737579167948677265543698044696895941600000000")*x_1^5 - ZZ("2212962118409540777636139481653983896606082752210036139528500000")*x_1^4 + ZZ("83343427603195335274689050203371826744335631347033627295380000")*x_1^3 - ZZ("4455083538531173585653387608123878667714263004789886085976728425")*x_1^2 + ZZ("8790723701979867940617716352259092695283280074229294598146348080")*x_1 + ZZ("13889925343373196298466410570440973429454153921683772288878965000")*x_3 - ZZ("1667930886634559159052902687235937686089504402779881982744812432"),
    109520000000000000*x_1^8 + 525475550000000000*x_1^7 + 783698651500000000*x_1^6 + 790767242590000000*x_1^5 - 10476732041700000*x_1^4 + 169669468523729500*x_1^3 - 16243383200380425*x_1^2 + 53082505589453640*x_1 - 18652836321211824
]

# Reducing the rational system is far cheaper than building the same
# polynomials again over the prime field.
prime = Int32(AlgebraicSolving.Nemo.rand_bits_prime(ZZ, 31))
Rp, _ = polynomial_ring(GF(prime), ["x_1", "x_2", "x_3"])
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
        println(io, join((NAME, 3, 3,
                          round(time_qq, digits = 6),
                          round(time_gf, digits = 6),
                          round(ratio, digits = 4)), ","))
    end
end
println("appended to $(csv)")
