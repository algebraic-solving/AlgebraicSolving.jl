# Benchmark for `equidimensional_decomposition` on the "jorge3" system.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AlgebraicSolving
using FileWatching.Pidfile: mkpidlock

const NAME = "jorge3"

R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])

F = [
    29250*x^3*y^2 - 113050*x^3*y*z - 217800*x^3*z^2 + 85800*x^2*y^3 - 195265*x^2*y^2*z - 516080*x^2*y^2 + 501390*x^2*y*z^2 - 797290*x^2*y*z - 119995*x^2*z^3 + 203790*x^2*z^2 + 81120*x*y^4 - 252285*x*y^3*z - 511894*x*y^3 + 392990*x*y^2*z^2 + 524534*x*y^2*z - 290756*x*y^2 + 168915*x*y*z^3 - 1072472*x*y*z^2 + 1353444*x*y*z - 16930*x*z^4 + 54320*x*z^3 - 40920*x*z^2 + 24960*y^5 - 134740*y^4*z - 79224*y^4 - 46540*y^3*z^2 + 644008*y^3*z + 53280*y^3 + 68560*y^2*z^3 - 151216*y^2*z^2 - 545568*y^2*z + 1054848*y^2 + 6300*y*z^4 - 148588*y*z^3 + 497312*y*z^2 - 450672*y*z,
    245440*x^3*y^2 - 218440*x^3*y*z + 118920*x^3*z^2 + 388362*x^2*y^3 - 1589246*x^2*y^2*z - 2531692*x^2*y^2 + 4270*x^2*y*z^2 + 2413552*x^2*y*z - 15002*x^2*z^3 - 275676*x^2*z^2 + 182858*x*y^4 - 767384*x*y^3*z - 2290656*x*y^3 - 1087960*x*y^2*z^2 + 2964612*x*y^2*z + 5231128*x*y^2 + 126912*x*y*z^3 + 875096*x*y*z^2 - 4239208*x*y*z - 17690*x*z^4 + 35908*x*z^3 + 167904*x*z^2 + 23400*y^5 + 317110*y^4*z - 263724*y^4 + 162978*y^3*z^2 - 62496*y^3*z + 2063160*y^3 - 62854*y^2*z^3 + 692408*y^2*z^2 - 820456*y^2*z - 3013536*y^2 + 2790*y*z^4 - 92240*y*z^3 - 668000*y*z^2 + 1952976*y*z - 1248*z^5 + 23940*z^4 - 14304*z^3 - 57168*z^2,
    51750*x^3*y^2 - 249090*x^3*y*z + 299460*x^3*z^2 + 87810*x^2*y^3 - 622808*x^2*y^2*z - 782298*x^2*y^2 + 1256768*x^2*y*z^2 + 3033480*x^2*y*z + 278606*x^2*z^3 - 513192*x^2*z^2 + 46578*x*y^4 - 283939*x*y^3*z - 773088*x*y^3 - 122385*x*y^2*z^2 + 3006255*x*y^2*z + 3363228*x*y^2 + 394487*x*y*z^3 - 862598*x*y*z^2 - 4510470*x*y*z + 100795*x*z^4 - 253931*x*z^3 + 236202*x*z^2 + 7560*y^5 + 62439*y^4*z - 154854*y^4 - 69311*y^3*z^2 - 106929*y^3*z + 1113174*y^3 - 105573*y^2*z^3 + 3408*y^2*z^2 - 1899918*y^2*z - 2934252*y^2 + 1951*y*z^4 - 264695*y*z^3 - 148908*y*z^2 + 1551420*y*z + 12774*z^5 - 37758*z^4 + 45444*z^3 - 42048*z^2,
    6624800*x^4*y + 4456800*x^4*z + 9428380*x^3*y^2 - 32123880*x^3*y*z - 17855680*x^3*y + 13889140*x^3*z^2 - 7898280*x^3*z - 4303208*x^2*y^3 - 45476486*x^2*y^2*z + 46018856*x^2*y^2 - 32385320*x^2*y*z^2 + 130710864*x^2*y*z - 116812448*x^2*y + 3087198*x^2*z^3 - 14697912*x^2*z^2 + 16999032*x^2*z - 11420032*x*y^4 + 17830726*x*y^3*z + 68337620*x*y^3 - 31601150*x*y^2*z^2 + 121570128*x*y^2*z - 168966056*x*y^2 - 2195978*x*y*z^3 + 60966800*x*y*z^2 - 256357336*x*y*z + 286261696*x*y - 294190*x*z^4 + 1627348*x*z^3 + 3887048*x*z^2 - 11929968*x*z - 4268160*y^5 + 25096600*y^4*z + 16075920*y^4 + 17577878*y^3*z^2 - 50526584*y^3*z - 329064*y^3 - 3261404*y^2*z^3 + 12307644*y^2*z^2 - 38129632*y^2*z + 53058480*y^2 - 77010*y*z^4 + 295640*y*z^3 - 30943256*y*z^2 + 133675520*y*z - 144710976*y - 51168*z^5 + 1083876*z^4 - 2549544*z^3 - 1170960*z^2 + 4687776*z,
    1035000*x^4*y - 2566800*x^4*z + 1092450*x^3*y^2 - 9976750*x^3*y*z - 1206210*x^3*y + 26275730*x^3*z^2 + 11993940*x^3*z - 1571280*x^2*y^3 - 3648703*x^2*y^2*z + 17592132*x^2*y^2 + 14898478*x^2*y*z^2 + 5834381*x^2*y*z - 60994644*x^2*y + 19447161*x^2*z^3 - 39250039*x^2*z^2 + 9592434*x^2*z - 2491872*x*y^4 + 10455011*x*y^3*z + 25674570*x*y^3 - 24472395*x*y^2*z^2 - 4269971*x*y^2*z - 124995966*x*y^2 + 6015037*x*y*z^3 - 18022208*x*y*z^2 - 44707334*x*y*z + 141802404*x*y + 5255985*x*z^4 - 19874311*x*z^3 + 25116044*x*z^2 - 12782724*x*z - 852480*y^5 + 6642444*y^4*z + 8663112*y^4 - 4430001*y^3*z^2 - 45928152*y^3*z - 16041132*y^3 - 9241278*y^2*z^3 + 2752721*y^2*z^2 + 9178620*y^2*z + 55929780*y^2 - 20959*y*z^4 - 2587003*y*z^3 + 2996442*y*z^2 + 33633948*y*z - 58222296*y + 523734*z^5 - 2595546*z^4 + 4959360*z^3 - 5450376*z^2 + 3447936*z,
    3036000*x^4*y - 7529280*x^4*z + 4754460*x^3*y^2 - 21886456*x^3*y*z - 55501776*x^3*y - 6125988*x^3*z^2 + 24113616*x^3*z + 1083858*x^2*y^3 + 2882442*x^2*y^2*z - 55614846*x^2*y^2 - 16238698*x^2*y*z^2 + 94986016*x^2*y*z + 305848404*x^2*y - 11096786*x^2*z^3 - 8598546*x^2*z^2 - 34761444*x^2*z - 1238286*x*y^4 + 11225644*x*y^3*z + 6469710*x*y^3 + 1209768*x*y^2*z^2 + 1278090*x*y^2*z + 138271764*x*y^2 + 11358052*x*y*z^3 + 83999142*x*y*z^2 - 78335080*x*y*z - 498731712*x*y - 4233002*x*z^4 + 5724082*x*z^3 + 16483548*x*z^2 + 19626192*x*z - 493560*y^5 - 1136814*y^4*z + 12653748*y^4 - 1771102*y^3*z^2 - 34027290*y^3*z - 71935092*y^3 + 448850*y^2*z^3 - 6002820*y^2*z^2 - 256656*y^2*z - 38576448*y^2 + 2065542*y*z^4 - 3495746*y*z^3 - 52845548*y*z^2 + 18346248*y*z + 236207952*y - 408084*z^5 + 2306412*z^4 + 2134032*z^3 - 6573024*z^2 - 7312032*z
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
        println(io, join((NAME, 3, 6,
                          round(time_qq, digits = 6),
                          round(time_gf, digits = 6),
                          round(ratio, digits = 4)), ","))
    end
end
println("appended to $(csv)")
