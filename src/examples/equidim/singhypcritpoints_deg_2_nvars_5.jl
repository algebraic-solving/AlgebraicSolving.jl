# Benchmark for `equidimensional_decomposition` on the "singhypcritpoints_deg_2_nvars_5" system.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AlgebraicSolving
using FileWatching.Pidfile: mkpidlock

const NAME = "singhypcritpoints_deg_2_nvars_5"

R, (x1, x2, x3, x4, x5) = polynomial_ring(QQ, ["x1", "x2", "x3", "x4", "x5"])

F = [
    -64547849*x1^3 + 394457886*x1^2*x2 - 116124080*x1^2*x3 - 52732783*x1^2*x4 + 133393756*x1^2*x5 - 39419961*x1^2 + 59372811*x1*x2^2 + 116659752*x1*x2*x3 + 319551738*x1*x2*x4 - 463702488*x1*x2*x5 + 376775106*x1*x2 + 54673956*x1*x3^2 + 117871213*x1*x3*x4 + 145316329*x1*x3*x5 + 148797844*x1*x3 + 35641532*x1*x4^2 - 18071759*x1*x4*x5 - 56930196*x1*x4 - 349538840*x1*x5^2 - 181415300*x1*x5 + 2766240*x1 + 122234976*x2^3 - 156183156*x2^2*x3 + 45507429*x2^2*x4 - 111721572*x2^2*x5 + 230330583*x2^2 + 126497318*x2*x3^2 + 176060180*x2*x3*x4 - 83846058*x2*x3*x5 - 37107472*x2*x3 + 205770902*x2*x4^2 - 516654600*x2*x4*x5 + 184030930*x2*x4 + 119684064*x2*x5^2 - 380465082*x2*x5 - 38726872*x2 - 15754388*x3^3 + 30567354*x3^2*x4 + 33530998*x3^2*x5 + 76315368*x3^2 + 54395307*x3*x4^2 + 26629291*x3*x4*x5 + 209095856*x3*x4 - 252048394*x3*x5^2 + 10941108*x3*x5 + 143985194*x3 + 48980115*x4^3 - 150140197*x4^2*x5 - 13581082*x4^2 - 89002780*x4*x5^2 - 94382932*x4*x5 - 114826269*x4 + 341246211*x5^3 - 107335469*x5^2 - 82074685*x5 - 56732004,
    68607752*x1^3 - 116124080*x1^2*x2 + 2676260*x1^2*x3 - 22485715*x1^2*x4 + 21016654*x1^2*x5 + 101626296*x1^2 + 58329876*x1*x2^2 + 109347912*x1*x2*x3 + 117871213*x1*x2*x4 + 145316329*x1*x2*x5 + 148797844*x1*x2 - 133896864*x1*x3^2 - 258804812*x1*x3*x4 + 174229796*x1*x3*x5 - 156250936*x1*x3 - 82728781*x1*x4^2 + 203619650*x1*x4*x5 - 148088774*x1*x4 + 29870113*x1*x5^2 + 195584094*x1*x5 + 101356988*x1 - 52061052*x2^3 + 126497318*x2^2*x3 + 88030090*x2^2*x4 - 41923029*x2^2*x5 - 18553736*x2^2 - 47263164*x2*x3^2 + 61134708*x2*x3*x4 + 67061996*x2*x3*x5 + 152630736*x2*x3 + 54395307*x2*x4^2 + 26629291*x2*x4*x5 + 209095856*x2*x4 - 252048394*x2*x5^2 + 10941108*x2*x5 + 143985194*x2 + 10456496*x3^3 - 32800932*x3^2*x4 + 26070960*x3^2*x5 + 23396904*x3^2 + 25489078*x3*x4^2 + 206319744*x3*x4*x5 - 125832448*x3*x4 + 116904382*x3*x5^2 - 4861520*x3*x5 + 83159616*x3 + 31790836*x4^3 + 40534157*x4^2*x5 - 55061438*x4^2 - 252519975*x4*x5^2 - 89679218*x4*x5 - 103814996*x4 - 261404167*x5^3 - 125501428*x5^2 - 238043640*x5 + 79558668,
    11474651*x1^3 - 52732783*x1^2*x2 - 22485715*x1^2*x3 + 26621818*x1^2*x4 - 9336737*x1^2*x5 - 35738237*x1^2 + 159775869*x1*x2^2 + 117871213*x1*x2*x3 + 71283064*x1*x2*x4 - 18071759*x1*x2*x5 - 56930196*x1*x2 - 129402406*x1*x3^2 - 165457562*x1*x3*x4 + 203619650*x1*x3*x5 - 148088774*x1*x3 - 73997829*x1*x4^2 + 217161602*x1*x4*x5 + 43742180*x1*x4 - 67784013*x1*x5^2 + 72770972*x1*x5 - 34164211*x1 + 15169143*x2^3 + 88030090*x2^2*x3 + 205770902*x2^2*x4 - 258327300*x2^2*x5 + 92015465*x2^2 + 30567354*x2*x3^2 + 108790614*x2*x3*x4 + 26629291*x2*x3*x5 + 209095856*x2*x3 + 146940345*x2*x4^2 - 300280394*x2*x4*x5 - 27162164*x2*x4 - 89002780*x2*x5^2 - 94382932*x2*x5 - 114826269*x2 - 10933644*x3^3 + 25489078*x3^2*x4 + 103159872*x3^2*x5 - 62916224*x3^2 + 95372508*x3*x4^2 + 81068314*x3*x4*x5 - 110122876*x3*x4 - 252519975*x3*x5^2 - 89679218*x3*x5 - 103814996*x3 + 62747716*x4^3 - 127601355*x4^2*x5 - 137102718*x4^2 - 183868948*x4*x5^2 + 78667960*x4*x5 + 12024124*x4 + 170886963*x5^3 + 197034765*x5^2 + 94025143*x5 + 63899514,
    -55970964*x1^3 + 133393756*x1^2*x2 + 21016654*x1^2*x3 - 9336737*x1^2*x4 - 79585218*x1^2*x5 - 102340163*x1^2 - 231851244*x1*x2^2 + 145316329*x1*x2*x3 - 18071759*x1*x2*x4 - 699077680*x1*x2*x5 - 181415300*x1*x2 + 87114898*x1*x3^2 + 203619650*x1*x3*x4 + 59740226*x1*x3*x5 + 195584094*x1*x3 + 108580801*x1*x4^2 - 135568026*x1*x4*x5 + 72770972*x1*x4 - 318855216*x1*x5^2 - 377338982*x1*x5 - 174253537*x1 - 37240524*x2^3 - 41923029*x2^2*x3 - 258327300*x2^2*x4 + 119684064*x2^2*x5 - 190232541*x2^2 + 33530998*x2*x3^2 + 26629291*x2*x3*x4 - 504096788*x2*x3*x5 + 10941108*x2*x3 - 150140197*x2*x4^2 - 178005560*x2*x4*x5 - 94382932*x2*x4 + 1023738633*x2*x5^2 - 214670938*x2*x5 - 82074685*x2 + 8690320*x3^3 + 103159872*x3^2*x4 + 116904382*x3^2*x5 - 2430760*x3^2 + 40534157*x3*x4^2 - 505039950*x3*x4*x5 - 89679218*x3*x4 - 784212501*x3*x5^2 - 251002856*x3*x5 - 238043640*x3 - 42533785*x4^3 - 183868948*x4^2*x5 + 39333980*x4^2 + 512660889*x4*x5^2 + 394069530*x4*x5 + 94025143*x4 + 857291136*x5^3 + 409516389*x5^2 + 419722578*x5 - 97582338,
    47686601*x1^4 - 64547849*x1^3*x2 + 68607752*x1^3*x3 + 11474651*x1^3*x4 - 55970964*x1^3*x5 + 145939547*x1^3 + 197228943*x1^2*x2^2 - 116124080*x1^2*x2*x3 - 52732783*x1^2*x2*x4 + 133393756*x1^2*x2*x5 - 39419961*x1^2*x2 + 1338130*x1^2*x3^2 - 22485715*x1^2*x3*x4 + 21016654*x1^2*x3*x5 + 101626296*x1^2*x3 + 13310909*x1^2*x4^2 - 9336737*x1^2*x4*x5 - 35738237*x1^2*x4 - 39792609*x1^2*x5^2 - 102340163*x1^2*x5 + 152510360*x1^2 + 19790937*x1*x2^3 + 58329876*x1*x2^2*x3 + 159775869*x1*x2^2*x4 - 231851244*x1*x2^2*x5 + 188387553*x1*x2^2 + 54673956*x1*x2*x3^2 + 117871213*x1*x2*x3*x4 + 145316329*x1*x2*x3*x5 + 148797844*x1*x2*x3 + 35641532*x1*x2*x4^2 - 18071759*x1*x2*x4*x5 - 56930196*x1*x2*x4 - 349538840*x1*x2*x5^2 - 181415300*x1*x2*x5 + 2766240*x1*x2 - 44632288*x1*x3^3 - 129402406*x1*x3^2*x4 + 87114898*x1*x3^2*x5 - 78125468*x1*x3^2 - 82728781*x1*x3*x4^2 + 203619650*x1*x3*x4*x5 - 148088774*x1*x3*x4 + 29870113*x1*x3*x5^2 + 195584094*x1*x3*x5 + 101356988*x1*x3 - 24665943*x1*x4^3 + 108580801*x1*x4^2*x5 + 21871090*x1*x4^2 - 67784013*x1*x4*x5^2 + 72770972*x1*x4*x5 - 34164211*x1*x4 - 106285072*x1*x5^3 - 188669491*x1*x5^2 - 174253537*x1*x5 + 101090304*x1 + 30558744*x2^4 - 52061052*x2^3*x3 + 15169143*x2^3*x4 - 37240524*x2^3*x5 + 76776861*x2^3 + 63248659*x2^2*x3^2 + 88030090*x2^2*x3*x4 - 41923029*x2^2*x3*x5 - 18553736*x2^2*x3 + 102885451*x2^2*x4^2 - 258327300*x2^2*x4*x5 + 92015465*x2^2*x4 + 59842032*x2^2*x5^2 - 190232541*x2^2*x5 - 19363436*x2^2 - 15754388*x2*x3^3 + 30567354*x2*x3^2*x4 + 33530998*x2*x3^2*x5 + 76315368*x2*x3^2 + 54395307*x2*x3*x4^2 + 26629291*x2*x3*x4*x5 + 209095856*x2*x3*x4 - 252048394*x2*x3*x5^2 + 10941108*x2*x3*x5 + 143985194*x2*x3 + 48980115*x2*x4^3 - 150140197*x2*x4^2*x5 - 13581082*x2*x4^2 - 89002780*x2*x4*x5^2 - 94382932*x2*x4*x5 - 114826269*x2*x4 + 341246211*x2*x5^3 - 107335469*x2*x5^2 - 82074685*x2*x5 - 56732004*x2 + 2614124*x3^4 - 10933644*x3^3*x4 + 8690320*x3^3*x5 + 7798968*x3^3 + 12744539*x3^2*x4^2 + 103159872*x3^2*x4*x5 - 62916224*x3^2*x4 + 58452191*x3^2*x5^2 - 2430760*x3^2*x5 + 41579808*x3^2 + 31790836*x3*x4^3 + 40534157*x3*x4^2*x5 - 55061438*x3*x4^2 - 252519975*x3*x4*x5^2 - 89679218*x3*x4*x5 - 103814996*x3*x4 - 261404167*x3*x5^3 - 125501428*x3*x5^2 - 238043640*x3*x5 + 79558668*x3 + 15686929*x4^4 - 42533785*x4^3*x5 - 45700906*x4^3 - 91934474*x4^2*x5^2 + 39333980*x4^2*x5 + 6012062*x4^2 + 170886963*x4*x5^3 + 197034765*x4*x5^2 + 94025143*x4*x5 + 63899514*x4 + 214322784*x5^4 + 136505463*x5^3 + 209861289*x5^2 - 97582338*x5 + 72716049
]

# Reducing the rational system is far cheaper than building the same
# polynomials again over the prime field.
prime = Int32(AlgebraicSolving.Nemo.rand_bits_prime(ZZ, 31))
Rp, _ = polynomial_ring(GF(prime), ["x1", "x2", "x3", "x4", "x5"])
Fp = [AlgebraicSolving.reduce_mod_p(f, Rp) for f in F]

# Decompose a tiny system first so that the timings below measure the
# computation rather than compilation. The solver specialises on the number of
# variables, so the warm up uses as many of them as the benchmark itself.
function warmup(K)
    _, w = polynomial_ring(K, ["w$(i)" for i in 1:5])
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
        println(io, join((NAME, 5, 5,
                          round(time_qq, digits = 6),
                          round(time_gf, digits = 6),
                          round(ratio, digits = 4)), ","))
    end
end
println("appended to $(csv)")
