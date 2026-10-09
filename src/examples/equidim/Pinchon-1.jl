# Benchmark for `equidimensional_decomposition` on the "Pinchon-1" system.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AlgebraicSolving
using FileWatching.Pidfile: mkpidlock

const NAME = "Pinchon-1"

R, (m16, m15, m14, m13, m12, m11, m10, m9, m8, m7, m6, m5, m4, m3, m2, m1, m0, h0, h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, A) = polynomial_ring(QQ, ["m16", "m15", "m14", "m13", "m12", "m11", "m10", "m9", "m8", "m7", "m6", "m5", "m4", "m3", "m2", "m1", "m0", "h0", "h1", "h2", "h3", "h4", "h5", "h6", "h7", "h8", "h9", "h10", "A"])

F = [
    m16 - h1 - 65536*h2 - 43046721*h3 - 4294967296*h4 - 152587890625*h5 - 2821109907456*h6 - 33232930569601*h7 - 281474976710656*h8 - 1853020188851841*h9 - 10000000000000000*h10,
    m15 - h1 - 32768*h2 - 14348907*h3 - 1073741824*h4 - 30517578125*h5 - 470184984576*h6 - 4747561509943*h7 - 35184372088832*h8 - 205891132094649*h9 - 1000000000000000*h10,
    m14 - h1 - 16384*h2 - 4782969*h3 - 268435456*h4 - 6103515625*h5 - 78364164096*h6 - 678223072849*h7 - 4398046511104*h8 - 22876792454961*h9 - 100000000000000*h10,
    m13 - h1 - 8192*h2 - 1594323*h3 - 67108864*h4 - 1220703125*h5 - 13060694016*h6 - 96889010407*h7 - 549755813888*h8 - 2541865828329*h9 - 10000000000000*h10,
    m12 - h1 - 4096*h2 - 531441*h3 - 16777216*h4 - 244140625*h5 - 2176782336*h6 - 13841287201*h7 - 68719476736*h8 - 282429536481*h9 - 1000000000000*h10,
    m11 - h1 - 2048*h2 - 177147*h3 - 4194304*h4 - 48828125*h5 - 362797056*h6 - 1977326743*h7 - 8589934592*h8 - 31381059609*h9 - 100000000000*h10,
    m10 - h1 - 1024*h2 - 59049*h3 - 1048576*h4 - 9765625*h5 - 60466176*h6 - 282475249*h7 - 1073741824*h8 - 3486784401*h9 - 10000000000*h10,
    m9 - h1 - 512*h2 - 19683*h3 - 262144*h4 - 1953125*h5 - 10077696*h6 - 40353607*h7 - 134217728*h8 - 387420489*h9 - 1000000000*h10,
    m8 - h1 - 256*h2 - 6561*h3 - 65536*h4 - 390625*h5 - 1679616*h6 - 5764801*h7 - 16777216*h8 - 43046721*h9 - 100000000*h10,
    m7 - h1 - 128*h2 - 2187*h3 - 16384*h4 - 78125*h5 - 279936*h6 - 823543*h7 - 2097152*h8 - 4782969*h9 - 10000000*h10,
    m6 - h1 - 64*h2 - 729*h3 - 4096*h4 - 15625*h5 - 46656*h6 - 117649*h7 - 262144*h8 - 531441*h9 - 1000000*h10,
    m5 - h1 - 32*h2 - 243*h3 - 1024*h4 - 3125*h5 - 7776*h6 - 16807*h7 - 32768*h8 - 59049*h9 - 100000*h10,
    m4 - h1 - 16*h2 - 81*h3 - 256*h4 - 625*h5 - 1296*h6 - 2401*h7 - 4096*h8 - 6561*h9 - 10000*h10,
    m3 - h1 - 8*h2 - 27*h3 - 64*h4 - 125*h5 - 216*h6 - 343*h7 - 512*h8 - 729*h9 - 1000*h10,
    m2 - h1 - 4*h2 - 9*h3 - 16*h4 - 25*h5 - 36*h6 - 49*h7 - 64*h8 - 81*h9 - 100*h10,
    m1 - h1 - 2*h2 - 3*h3 - 4*h4 - 5*h5 - 6*h6 - 7*h7 - 8*h8 - 9*h9 - 10*h10,
    m0 - h0 - h1 - h2 - h3 - h4 - h5 - h6 - h7 - h8 - h9 - h10,
    h0 - h1 + h2 - h3 + h4 - h5 + h6 - h7 + h8 - h9 + h10,
    m1 - A,
    m0 - 1,
    m2 - A^2,
    m3 - A^3,
    m4 - A^4,
    m6 - 6*m5*A + 5*A^6,
    m8 - 8*m7*A + 112*m5*A^3 - 105*A^8,
    m10 - 10*m9*A + 240*m7*A^3 - 126*m5^2 - 3780*m5*A^5 + 3675*A^10,
    m12 - 12*m11*A + 440*m9*A^3 - 792*m7*m5 - 11880*m7*A^5 + 24948*m5^2*A^2 + 166320*m5*A^7 - 179025*A^12,
    m14 - 14*m13*A + 728*m11*A^3 - 2002*m9*m5 - 30030*m9*A^5 - 1716*m7^2 + 216216*m7*m5*A^2 + 720720*m7*A^7 - 4162158*m5^2*A^4 - 7777770*m5*A^9 + 11036025*A^14,
    m16 - 16*m15*A + 1120*m13*A^3 - 4368*m11*m5 - 65520*m11*A^5 - 11440*m9*m7 + 720720*m9*m5*A^2 + 2402400*m9*A^7 + 617760*m7^2*A^2 - 47567520*m7*m5*A^4 - 44444400*m7*A^9 + 6054048*m5^3*A + 711350640*m5^2*A^6 + 151351200*m5*A^11 - 780404625*A^16
]

# Reducing the rational system is far cheaper than building the same
# polynomials again over the prime field.
prime = Int32(AlgebraicSolving.Nemo.rand_bits_prime(ZZ, 31))
Rp, _ = polynomial_ring(GF(prime), ["m16", "m15", "m14", "m13", "m12", "m11", "m10", "m9", "m8", "m7", "m6", "m5", "m4", "m3", "m2", "m1", "m0", "h0", "h1", "h2", "h3", "h4", "h5", "h6", "h7", "h8", "h9", "h10", "A"])
Fp = [AlgebraicSolving.reduce_mod_p(f, Rp) for f in F]

# Decompose a tiny system first so that the timings below measure the
# computation rather than compilation. The solver specialises on the number of
# variables, so the warm up uses as many of them as the benchmark itself.
function warmup(K)
    _, w = polynomial_ring(K, ["w$(i)" for i in 1:29])
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
        println(io, join((NAME, 29, 29,
                          round(time_qq, digits = 6),
                          round(time_gf, digits = 6),
                          round(ratio, digits = 4)), ","))
    end
end
println("appended to $(csv)")
