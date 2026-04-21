import msolve_jll: libmsolve

export real_solutions, rational_solutions, rational_parametrization

@doc Markdown.doc"""
    _get_rational_parametrization(nr::Int32, lens::Vector{Int32}, cfs::Ptr{BigInt})

Construct the rational parametrization of the solution set computed via msolve.

**Note**: This is an internal function and should only be used inside `msolve()`.
"""
function _get_rational_parametrization(
        nr::Int32,
        lens::Vector{Int32},
        cfs::Vector{BigInt},
        cfs_lf::Vector{BigInt},
        nr_vars::Int32
    )
    if !isempty(cfs_lf)
        lf = [ZZRingElem(cfs_lf[i]) for i in 1:nr_vars]
    else
        lf = ZZRingElem[]
    end
    C, x  = polynomial_ring(QQ,"x")
    ctr   = 0

    elim  = C([cfs[i] for i in 1:lens[1]])
    ctr += lens[1]

    denom = C([cfs[i+ctr] for i in 1:lens[2]])
    ctr +=  lens[2]

    size  = nr-2
    p = Vector{QQPolyRingElem}(undef, size)
    k = 1
    for i in 3:nr
        p[k]  = C([cfs[j+ctr] for j in 1:lens[i]-1])
        # multiply parametrization polynomial directly with
        # corresponding coefficients
        p[k]  *=   (-1) // ZZRingElem(cfs[lens[i]+ctr])
        ctr   +=  lens[i]
        k     +=  1
    end

    return lf, elim, denom, p
end

@doc Markdown.doc"""
    _core_msolve_array(lens::Vector{Int32}, cfs::Union{Vector{Int32},Vector{BigInt}}, exps::Vector{Int32}, variable_names::Vector{String}, field_char::Int; <keyword arguments>)

Compute a rational parametrization and the real solutions of the ideal via msolve.

**Note**: Both the input and output of this function are given as flattened arrays of integers and strings instead of polynomials.

# Arguments
- `lens::Vector{Int32}`: vector of number of terms per generator.
- `cfs::Union{Vector{Int32},Vector{BigInt}}`: vector of coefficients of all generators, concatenated. When working over the rationals, the coefficients are given as pairs of numerators and denominators, i.e. `cfs[2*i-1]` is the numerator and `cfs[2*i]` the denominator of the `i`-th coefficient.
- `exps::Vector{Int32}`: vector of exponent vectors of all generators, concatenated.
- `variable_names::Vector{String}`: vector of variable names.
- `field_char::Int`: characteristic of the ground field, `0` for the rationals.
- `initial_hts::Int=17`: initial hash table size `log_2`.
- `nr_thrds::Int=1`: number of threads for parallel linear algebra.
- `max_nr_pairs::Int=0`: maximal number of pairs per matrix, only bounded by minimal degree if `0`.
- `la_option::Int=2`: linear algebra option: exact sparse-dense (`1`), exact sparse (`2`, default), probabilistic sparse-dense (`42`), probabilistic sparse(`44`).
- `get_param::Int=1`: parametrization option: no (`0`), yes (`1`, default), only rational parametrization, no solutions are computed (`2`).
- `info_level::Int=0`: info level printout: off (`0`, default), summary (`1`), detailed (`2`).
- `precision::Int=32`: bit precision for the computed solutions.
```

**Note**: This is an internal function.
"""
function _core_msolve_array(
    lens::Vector{Int32},
    cfs::Union{Vector{Int32},Vector{BigInt}},
    exps::Vector{Int32},
    variable_names::Vector{String},
    field_char::Int;
    initial_hts::Int=17,
    nr_thrds::Int=1,
    max_nr_pairs::Int=0,
    la_option::Int=2,
    get_param::Int=1,
    info_level::Int=0,
    precision::Int=32
)::Tuple{Vector{Int32},Vector{String},Vector{BigInt},Vector{BigInt},Vector{BigInt},Vector{Int32}}
    if field_char != 0
        error("At the moment we only support the rationals as ground field.")
    end

    nr_gens = length(lens)
    nr_vars = length(exps) ÷ sum(lens)

    #= add new variables and linear forms, =#
    genericity_handling = 2

    reset_ht = 0
    print_gb = 0

    mon_order = 0
    elim_block_size = 0
    use_signatures = 0

    res_ld = Ref(Cint(0))
    res_nr_vars = Ref(Cint(0))
    res_dim = Ref(Cint(0))
    res_dquot = Ref(Cint(0))
    nb_sols = Ref(Cint(0))
    res_len = Ref(Ptr{Cint}(0))
    res_vnames = Ref(Ptr{Ptr{Cchar}}(0))
    res_cf_lf = Ref(Ptr{Cvoid}(0))
    res_cf = Ref(Ptr{Cvoid}(0))
    sols_num = Ref(Ptr{Cvoid}(0))
    sols_den = Ref(Ptr{Cint}(0))

    ccall((:msolve_julia, libmsolve), Cvoid,
        (Ptr{Nothing}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},
            Ptr{Ptr{Cint}}, Ptr{Ptr{Ptr{Cchar}}}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cint}, Ptr{Cvoid}, Ptr{Ptr{Cint}},
            Ptr{Cint}, Ptr{Cint}, Ptr{Cvoid}, Ptr{Ptr{Cchar}}, Ptr{Cchar}, Cint, Cint,
            Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint,
            Cint, Cint, Cint, Cint, Cint, Cint),
        cglobal(:jl_malloc), res_ld, res_nr_vars, res_dim, res_dquot,
        res_len, res_vnames, res_cf_lf, res_cf, nb_sols, sols_num, sols_den,
        lens, exps, cfs, variable_names, "/dev/null", field_char, mon_order,
        elim_block_size, nr_vars, nr_gens, initial_hts, nr_thrds, max_nr_pairs, reset_ht, la_option,
        use_signatures, print_gb, get_param, genericity_handling, precision, info_level)

    # convert to julia array, also give memory management to julia
    jl_ld = res_ld[]
    jl_rp_nr_vars = res_nr_vars[]
    jl_dim = res_dim[]
    jl_dquot = res_dquot[]

    if jl_dim > 0
        error("Dimension of ideal is greater than zero, no solutions provided.")
    end

    if jl_dquot == 0
        return Int32[], String[], BigInt[], BigInt[], BigInt[], Int32[]
    end

    # get rational parametrization
    jl_len = Vector{Int32}(Base.unsafe_wrap(Array, res_len[], jl_ld))

    jl_vnames_ptrs = Base.unsafe_wrap(Array, res_vnames[], jl_rp_nr_vars)
    jl_vnames = [unsafe_string(jl_vnames_ptrs[i]) for i in 1:jl_rp_nr_vars]

    jl_cf_lf_ptr = reinterpret(Ptr{BigInt}, res_cf_lf[])
    if jl_cf_lf_ptr != Ptr{BigInt}(0)
        jl_cf_lf = BigInt[deepcopy(unsafe_load(jl_cf_lf_ptr, i)) for i in 1:jl_rp_nr_vars]
    else
        jl_cf_lf = BigInt[]
    end
    jl_cf_ptr = reinterpret(Ptr{BigInt}, res_cf[])
    jl_cf = BigInt[deepcopy(unsafe_load(jl_cf_ptr, i)) for i in 1:sum(jl_len)]

    jl_nb_sols = nb_sols[]

    # get solutions
    len = 2 * jl_nb_sols * nr_vars
    jl_sols_num_ptr = reinterpret(Ptr{BigInt}, sols_num[])
    jl_sols_num = BigInt[deepcopy(unsafe_load(jl_sols_num_ptr, i)) for i in 1:len]
    jl_sols_den_ptr = reinterpret(Ptr{Int32}, sols_den[])
    jl_sols_den = Vector{Int32}(Base.unsafe_wrap(Array, jl_sols_den_ptr, len))

    ccall((:free_msolve_julia_result_data, libmsolve), Nothing,
        (Ptr{Nothing}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cvoid}},
            Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cint}}, Cint, Cint, Cint),
        cglobal(:jl_free), res_len, res_cf, sols_num, sols_den,
        jl_ld, jl_nb_sols, field_char)

    return jl_len, jl_vnames, jl_cf_lf, jl_cf, jl_sols_num, jl_sols_den
end

@doc Markdown.doc"""
    _core_msolve(I::Ideal{T} where T <: MPolyRingElem, <keyword arguments>)

Compute a rational parametrization and the real solutions of the ideal `I` via msolve.

**Note**: This is an internal function.
"""
function _core_msolve(
        I::Ideal{T} where T <: MPolyRingElem;     # input generators
        initial_hts::Int=17,                  # hash table size, default 2^17
        nr_thrds::Int=1,                      # number of threads
        max_nr_pairs::Int=0,                  # number of pairs maximally chosen
                                              # in symbolic preprocessing
        la_option::Int=2,                     # linear algebra option
        get_param::Int=1,                     # get rational parametrization
        info_level::Int=0,                    # info level for print outs
        precision::Int=32                     # precision of the solution set
        )

    F = I.gens
    if (is_zero(F))
        error("Dimension of ideal is greater than zero, no solutions provided.")
    end
    R = first(F).parent
    nr_vars     = nvars(R)
    nr_gens     = length(F)
    field_char  = Int(characteristic(R))

    variable_names = map(string, R.S)

    if field_char != 0
        error("At the moment we only support the rationals as ground field.")
    end
    # convert Singular ideal to flattened arrays of ints
    lens, cfs, exps, nr_gens = _convert_to_msolve(F)

    jl_len, jl_vnames, jl_cf_lf, jl_cf, jl_sols_num, jl_sols_den =
        _core_msolve_array(lens, cfs, exps, variable_names, field_char;
            initial_hts, nr_thrds, max_nr_pairs, la_option, get_param, info_level, precision)

    # convert to julia array, also give memory management to julia
    jl_ld = Int32(length(jl_len))
    jl_rp_nr_vars = Int32(length(jl_vnames))
    jl_nb_sols = length(jl_sols_num) ÷ (2 * nr_vars)

    I.dim = 0

    # get rational parametrization
    nterms  = 0

    if jl_ld == 0
        I.dim = -1
        C,x = polynomial_ring(QQ,"x")
        I.rat_param = RationalParametrization(Symbol[], ZZRingElem[], C(-1), C(-1), QQPolyRingElem[])
        I.real_sols = QQFieldElem[]
        I.inter_sols = Vector{QQFieldElem}[]
        return I.rat_param, I.real_sols
    end
    [nterms += jl_len[i] for i=1:jl_ld]

    vsymbols  = [Symbol(jl_vnames[i]) for i in 1:jl_rp_nr_vars]
    #= get possible variable permutation, ignoring additional variables=#
    perm      = indexin(R.S, vsymbols)

    rat_param = _get_rational_parametrization(jl_ld, jl_len,
                                              jl_cf, jl_cf_lf, jl_rp_nr_vars)


    I.rat_param = RationalParametrization(vsymbols, rat_param[1],rat_param[2],
                                          rat_param[3], rat_param[4])
    if get_param == 2
        return rat_param, undef
    end

    if jl_nb_sols == 0
        I.real_sols = QQFieldElem[]
        I.inter_sols = Vector{QQFieldElem}[]
        return rat_param, Vector{QQFieldElem}[]
    end

    #= solutions are returned as intervals, i.e. a minimum and a maximum entry for
     = the numerator and denominator; we also sum up and divide by 2 =#

    inter_solutions = Vector{Vector{Vector{QQFieldElem}}}(undef, jl_nb_sols)
    solutions = Vector{Vector{QQFieldElem}}(undef, jl_nb_sols)

    len = 2*jl_nb_sols*nr_vars
    i = 1
    k = 1
    while i <= len
        j = 1
        inter_tmp = Vector{Vector{QQFieldElem}}(undef, nr_vars)
        tmp = Vector{QQFieldElem}(undef, nr_vars)
        while j <= nr_vars
            inter_tmp_coeff = Vector{QQFieldElem}(undef, 2)
            inter_tmp_coeff[1] = QQFieldElem(jl_sols_num[i]) >> Int64(jl_sols_den[i])
            inter_tmp_coeff[2] = QQFieldElem(jl_sols_num[i+1]) >> Int64(jl_sols_den[i+1])
            inter_tmp[perm[j]] = inter_tmp_coeff
            tmp[perm[j]] = sum(inter_tmp_coeff) >> 1
            i += 2
            j += 1
        end
        inter_solutions[k] = inter_tmp
        solutions[k] = tmp
        k += 1
    end
    I.inter_sols = inter_solutions
    I.real_sols = solutions

    return rat_param, solutions
end

@doc Markdown.doc"""
    rational_parametrization(I::Ideal{T} where T <: MPolyRingElem, <keyword arguments>)

Given an ideal `I` with a finite solution set over the complex numbers, return
the rational parametrization of the ideal with a given precision (default 32 bits).

**Note**: At the moment only QQ is supported as ground field. If the dimension of the ideal
is greater than zero an ErrorException is thrown.

# Arguments
- `I::Ideal{T} where T <: MPolyRingElem`: input generators.
- `initial_hts::Int=17`: initial hash table size `log_2`.
- `nr_thrds::Int=1`: number of threads for parallel linear algebra.
- `max_nr_pairs::Int=0`: maximal number of pairs per matrix, only bounded by minimal degree if `0`.
- `la_option::Int=2`: linear algebra option: exact sparse-dense (`1`), exact sparse (`2`, default), probabilistic sparse-dense (`42`), probabilistic sparse(`44`).
- `info_level::Int=0`: info level printout: off (`0`, default), summary (`1`), detailed (`2`).
- `precision::Int=32`: bit precision for the computed solutions.

# Examples
```jldoctest
julia> using AlgebraicSolving

julia> R,(x1,x2,x3) = polynomial_ring(QQ, ["x1","x2","x3"])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x1, x2, x3])

julia> I = Ideal([x1+2*x2+2*x3-1, x1^2+2*x2^2+2*x3^2-x1, 2*x1*x2+2*x2*x3-x2])
QQMPolyRingElem[x1 + 2*x2 + 2*x3 - 1, x1^2 - x1 + 2*x2^2 + 2*x3^2, 2*x1*x2 + 2*x2*x3 - x2]

julia> rational_parametrization(I)
AlgebraicSolving.RationalParametrization([:x1, :x2, :x3], ZZRingElem[], 84*x^4 - 40*x^3 + x^2 + x, 336*x^3 - 120*x^2 + 2*x + 1, Nemo.QQPolyRingElem[184*x^3 - 80*x^2 + 4*x + 1, 36*x^3 - 18*x^2 + 2*x])
```
"""
function rational_parametrization(
        I::Ideal{T} where T <: MPolyRingElem;     # input generators
        initial_hts::Int=17,                  # hash table size, default 2^17
        nr_thrds::Int=1,                      # number of threads
        max_nr_pairs::Int=0,                  # number of pairs maximally chosen
                                              # in symbolic preprocessing
        la_option::Int=2,                     # linear algebra option
        info_level::Int=0,                    # info level for print outs
        precision::Int=32                     # precision of the solution set
        )

    isdefined(I, :rat_param) ||
    _core_msolve(I,
                 initial_hts = initial_hts,
                 nr_thrds = nr_thrds,
                 max_nr_pairs = max_nr_pairs,
                 la_option = la_option,
                 info_level = info_level,
                 get_param = 2, # only rational parametrization, no solutions are computed
                 precision = precision)

    return I.rat_param
end


@doc Markdown.doc"""
    rational_solutions(I::Ideal{T} where T <: MPolyRingElem, <keyword arguments>)

Given an ideal `I` with a finite solution set over the complex numbers, return
the rational roots of the ideal. 

# Arguments
- `I::Ideal{T} where T <: MPolyRingElem`: input generators.
- `initial_hts::Int=17`: initial hash table size `log_2`.
- `nr_thrds::Int=1`: number of threads for parallel linear algebra.
- `max_nr_pairs::Int=0`: maximal number of pairs per matrix, only bounded by minimal degree if `0`.
- `la_option::Int=2`: linear algebra option: exact sparse-dense (`1`), exact sparse (`2`, default), probabilistic sparse-dense (`42`), probabilistic sparse(`44`).
- `info_level::Int=0`: info level printout: off (`0`, default), summary (`1`), detailed (`2`).
- `precision::Int=32`: bit precision for the computed solutions.

# Examples
```jldoctest
julia> using AlgebraicSolving

julia> R,(x1,x2,x3) = polynomial_ring(QQ, ["x1","x2","x3"])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x1, x2, x3])

julia> I = Ideal([x1+2*x2+2*x3-1, x1^2+2*x2^2+2*x3^2-x1, 2*x1*x2+2*x2*x3-x2])
QQMPolyRingElem[x1 + 2*x2 + 2*x3 - 1, x1^2 - x1 + 2*x2^2 + 2*x3^2, 2*x1*x2 + 2*x2*x3 - x2]

julia> rat_sols = rational_solutions(I)
2-element Vector{Vector{QQFieldElem}}:
 [1//3, 0, 1//3]
 [1, 0, 0]

julia> map(r->map(p->evaluate(p, r), I.gens), rat_sols)
2-element Vector{Vector{QQFieldElem}}:
 [0, 0, 0]
 [0, 0, 0]
```
"""
function rational_solutions(
        I::Ideal{T} where T <: MPolyRingElem;     # input generators
        initial_hts::Int=17,                  # hash table size, default 2^17
        nr_thrds::Int=1,                      # number of threads
        max_nr_pairs::Int=0,                  # number of pairs maximally chosen
                                              # in symbolic preprocessing
        la_option::Int=2,                     # linear algebra option
        info_level::Int=0,                    # info level for print outs
        precision::Int=32                     # precision of the solution set
        )
    isdefined(I, :rat_param) ||
    _core_msolve(I,
                 initial_hts = initial_hts,
                 nr_thrds = nr_thrds,
                 max_nr_pairs = max_nr_pairs,
                 la_option = la_option,
                 info_level = info_level,
                 get_param = 2, # only rational parametrization, no solutions are computed
                 precision = precision)
    param_t = I.rat_param

    # ge(t variable permutation
    R = parent(first(I.gens))
    perm = indexin(R.S, param_t.vars)
    nvars = length(param_t.vars)
    lpol = filter(l->degree(l) == 1, first.(collect(factor(param_t.elim))))
    nb = length(lpol)

    rat_elim = [-coeff(l, 0)// coeff(l, 1) for l in lpol]
    rat_den = map(l->evaluate(param_t.denom, l), rat_elim)
    rat_num = map(r->map(l->evaluate(l, r), param_t.param), rat_elim)

    rat_sols = Vector{Vector{QQFieldElem}}(undef, nb)

    if length(param_t.vars) == parent(I).nvars

      for i in 1:nb
        rat_sols[i] = Vector{QQFieldElem}(undef, nvars)
        for j in 1:(nvars-1)
           rat_sols[i][j] = rat_num[i][j] // rat_den[i]
        end
        rat_sols[i][nvars] = rat_elim[i]
      end

    else

      for i in 1:nb
        rat_sols[i] = Vector{QQFieldElem}(undef, nvars - 1)
        for j in 1:(nvars-1)
           rat_sols[i][j] = rat_num[i][j] // rat_den[i]
        end
      end

    end

    I.rat_sols = [r[perm] for r in rat_sols]

    return I.rat_sols
end


@doc Markdown.doc"""
    real_solutions(I::Ideal{T} where T <: MPolyRingElem, <keyword arguments>)

Given an ideal `I` with a finite solution set over the complex numbers, return
the real roots of the ideal with a given precision (default 32 bits).

**Note**: At the moment only QQ is supported as ground field. If the dimension of the ideal
is greater than zero an ErrorException is thrown.

# Arguments
- `I::Ideal{T} where T <: MPolyRingElem`: input generators.
- `initial_hts::Int=17`: initial hash table size `log_2`.
- `nr_thrds::Int=1`: number of threads for parallel linear algebra.
- `max_nr_pairs::Int=0`: maximal number of pairs per matrix, only bounded by minimal degree if `0`.
- `la_option::Int=2`: linear algebra option: exact sparse-dense (`1`), exact sparse (`2`, default), probabilistic sparse-dense (`42`), probabilistic sparse(`44`).
- `info_level::Int=0`: info level printout: off (`0`, default), summary (`1`), detailed (`2`).
- `precision::Int=32`: bit precision for the computed solutions.
- `interval::Bool=false`: return real solutions as intervals.

# Examples
```jldoctest
julia> using AlgebraicSolving

julia> R,(x1,x2,x3) = polynomial_ring(QQ, ["x1","x2","x3"])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x1, x2, x3])

julia> I = Ideal([x1+2*x2+2*x3-1, x1^2+2*x2^2+2*x3^2-x1, 2*x1*x2+2*x2*x3-x2])
QQMPolyRingElem[x1 + 2*x2 + 2*x3 - 1, x1^2 - x1 + 2*x2^2 + 2*x3^2, 2*x1*x2 + 2*x2*x3 - x2]

julia> real_solutions(I)
4-element Vector{Vector{QQFieldElem}}:
 [1454069069139182487//2305843009213693952, 1454069069139182487//4611686018427387904, -154187553040555781639//1180591620717411303424]
 [1, 0, 0]
 [8357884734418025835//36893488147419103232, 4178942367209012917//36893488147419103232, 40355437357166103123//147573952589676412928]
 [6148914691236517205//18446744073709551616, 123527304065019319//85070591730234615865843651857942052864, 196765270119568550571//590295810358705651712]

julia> real_solutions(I, interval=true)
4-element Vector{Vector{Vector{QQFieldElem}}}:
 [[11632552553113459895//18446744073709551616, 11632552553113459897//18446744073709551616], [5816276276556729947//18446744073709551616, 5816276276556729949//18446744073709551616], [-19273444130069472705//147573952589676412928, -77093776520277890819//590295810358705651712]]
 [[1, 1], [0, 0], [0, 0]]
 [[2089471183604506445//9223372036854775808, 4178942367209012945//18446744073709551616], [1044735591802253225//9223372036854775808, 2089471183604506467//18446744073709551616], [20177718678583051561//73786976294838206464, 10088859339291525781//36893488147419103232]]
 [[6148914691236517197//18446744073709551616, 6148914691236517213//18446744073709551616], [-9388075108941468237//170141183460469231731687303715884105728, 9882184325201545513//170141183460469231731687303715884105728], [98382635059784275285//295147905179352825856, 49191317529892137643//147573952589676412928]]
```
"""
function real_solutions(
        I::Ideal{T} where T <: MPolyRingElem;     # input generators
        initial_hts::Int=17,                  # hash table size, default 2^17
        nr_thrds::Int=1,                      # number of threads
        max_nr_pairs::Int=0,                  # number of pairs maximally chosen
                                              # in symbolic preprocessing
        la_option::Int=2,                     # linear algebra option
        info_level::Int=0,                    # info level for print outs
        precision::Int=32,                    # precision of the solution set
        interval::Bool=false                  # return real solutions as intervals
        )

    isdefined(I, :real_sols) ||
    _core_msolve(I,
                 initial_hts = initial_hts,
                 nr_thrds = nr_thrds,
                 max_nr_pairs = max_nr_pairs,
                 la_option = la_option,
                 info_level = info_level,
                 precision = precision)

    if interval
        return I.inter_sols
    else
        return I.real_sols
    end
end
