@doc Markdown.doc"""
    rational_parametrization(I::ParametricIdeal{K}; <keyword arguments>) -> (Vector{Symbol}, Vector{MPoly{K}}, MPoly{K}, MPoly{K}, Vector{MPoly{K}})

Computes the parametric rational parametrization associated with the parametric ideal `I`. The rational parametrization is computed by evaluating the ideal at a number of parameter values and interpolating the coefficients of the resulting rational parametrizations.

# Arguments
- `retry::Int=10`: the maximum number of consecutive failures allowed when computing the rational parametrization or interpolating the coefficients.
- `nr_thrds::Int=1`: the number of threads to use for parallel computations.
- `worker_pool::AbstractWorkerPool=default_worker_pool()`: the worker pool to use for parallel computations.
- `show_progress::Bool=true`: whether to show a progress bar while computing the rational parametrization.
"""
function rational_parametrization(
    I::ParametricIdeal{K} where K<:FracFieldElem;
    retry::Int=10,
    nr_thrds::Int=1,
    worker_pool::AbstractWorkerPool=default_worker_pool(),
    show_progress::Bool=true
)
    n = I.num_vars
    R = parent(I.gens[1])
    R₂ = base_ring(I.base_frac_field)
    R₄, _ = polynomial_ring(fraction_field(R₂), symbols(R)[1:n], internal_ordering=:degrevlex)
    param_shift = nothing
    function bb(v)
        gb = groebner_basis(I, v; worker_pool=worker_pool)
        param = AlgebraicSolving.rational_parametrization(AlgebraicSolving.Ideal(gb))
        vars = param.vars
        cfs_lf = param.cfs_lf
        elim = param.elim / leading_coefficient(param.elim)
        denom = param.denom / leading_coefficient(param.denom)
        params = map(p -> p / leading_coefficient(param.denom), param.param)
        param = (vars, cfs_lf, elim, denom, params)
        if !isnothing(param_shift)
            if param[1] != param_shift[1] ||
               param[2] != param_shift[2] ||
               degree(param[3]) != degree(param_shift[3]) ||
               degree(param[4]) != degree(param_shift[4]) ||
               degree.(param[5]) != degree.(param_shift[5])
                error("Error evaluating black box function: the rational parametrization has different monomial structure at different parameter values.")
            end
        end
        param
    end
    param_shift = bb(QQ.(I.shift))
    param_vars = param_shift[1]
    param_cfs_lf = param_shift[2]
    if length(param_cfs_lf) > n
        error("Error: the last variable in the rational parametrization is not a known variable: is the ideal radical?")
    end
    param_var_index = findfirst(i -> symbols(R₄)[i] == param_vars[end], 1:n)
    param_var = gens(R₄)[param_var_index]
    param_elim = zero(R₄)
    for i in 0:degree(param_shift[3])
        c = Interpolation.cuyt_lee(
            R₂, v -> coeff(bb(v)[3], i);
            initial_shift=I.shift,
            retry=retry,
            nr_thrds=nr_thrds,
            show_progress=show_progress,
            desc="Parametric rational parametrization (eliminating polynomial), deg=($(degree(param_shift[3]))), index=($i)"
        )
        param_elim += R₄(c) * param_var^i
    end
    param_denom = zero(R₄)
    for i in 0:degree(param_shift[4])
        c = Interpolation.cuyt_lee(
            R₂, v -> coeff(bb(v)[4], i);
            initial_shift=I.shift,
            retry=retry,
            nr_thrds=nr_thrds,
            show_progress=show_progress,
            desc="Parametric rational parametrization (denominator), deg=($(degree(param_shift[4]))), index=($i)"
        )
        param_denom += R₄(c) * param_var^i
    end
    param_params = [zero(R₄) for _ in 1:(n-1)]
    len = n - 1
    degs = degree.(param_shift[5])
    for i in 1:len
        for j in 0:degs[i]
            c = Interpolation.cuyt_lee(
                R₂, v -> coeff(bb(v)[5][i], j);
                initial_shift=I.shift,
                retry=retry,
                nr_thrds=nr_thrds,
                show_progress=show_progress,
                desc="Parametric rational parametrization (parametrizations), size=($len,$(maximum(degs))), index=($i,$j)"
            )
            param_params[i] += R₄(c) * param_var^j
        end
    end
    (param_vars, param_cfs_lf, param_elim, param_denom, param_params)
end
