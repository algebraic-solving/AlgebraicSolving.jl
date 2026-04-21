module RAG

using Nemo

export has_real_solutions, points_per_components

function _flint_get_str_pretty(p::QQMPolyRingElem)::String
    R = parent(p)
    x = map(var -> Base.unsafe_convert(Cstring, string(var)), gens(R))
    str = @ccall Nemo.libflint.fmpq_mpoly_get_str_pretty(p::Ref{QQMPolyRingElem}, x::Ptr{Ref{Cstring}}, R::Ref{QQMPolyRing})::Cstring
    unsafe_string(str)
end

function _get_str(f::Vector{QQMPolyRingElem})::String
    join([_flint_get_str_pretty(p) for p in f], ", ")
end

function _call_procedure_with_maple(
    procedure::String;
    savelibname::String=joinpath(homedir(), "libs"),
    maple_path::String="maple",
    msolve_path::String="msolve",
)::String
    if !isfile(joinpath(savelibname, "MSolve.mla"))
        error("Could not find MSolve.mla in $savelibname. Please set savelibname to the correct path.")
    end
    if isnothing(Sys.which(maple_path))
        error("Could not find Maple executable. Please set maple_path to the correct path.")
    end
    if isnothing(Sys.which(msolve_path))
        error("Could not find MSolve executable. Please set msolve_path to the correct path.")
    end
    input = tempname()
    output = tempname()
    open(input, "w") do io
        println(io, "savelibname := \"$savelibname\":")
        println(io, "libname := savelibname, libname:")
        println(io, "MSolve:-SetMSolvePath(\"$msolve_path\"):")
        println(io, "f:=proc()")
        println(io, procedure)
        println(io, "end proc:")
        println(io, "fd := fopen(\"$output\", WRITE):")
        println(io, "fprintf(fd, \"%a\\n\", f()):")
        println(io, "fclose(fd):")
    end
    command = Cmd(`$maple_path -q $input`)
    run(command)
    line = open(output, "r") do io
        readline(io)
    end
    line
end

function has_real_solutions(
    eqs::Vector{QQMPolyRingElem},
    pos::Vector{QQMPolyRingElem},
    ineqs::Vector{QQMPolyRingElem};
    savelibname::String=joinpath(homedir(), "libs"),
    maple_path::String="maple",
    msolve_path::String="msolve",
    info_level::Int=0
)::Vector{Vector{Vector{QQFieldElem}}}
    procedure = """
        local vars, res:
        vars := [$(_get_str(gens(parent([eqs; pos; ineqs][1]))))]:
        res := RAG:-HasRealSolutions([$(_get_str(eqs))], [$(_get_str(pos))], [$(_get_str(ineqs))], {"verb"=$info_level}):
        return map(point -> subs(point, vars), res):
    """
    line = _call_procedure_with_maple(procedure; savelibname, maple_path, msolve_path)
    line = replace(line, r"(-?\d+)/(\d+)" => s"QQ(\1)/QQ(\2)")
    line = replace(line, r"(-?\d+)" => s"QQ(\1)")
    res = eval(Meta.parse(line))
    res
end

function points_per_components(
    eqs::Vector{QQMPolyRingElem},
    pos::Vector{QQMPolyRingElem},
    ineqs::Vector{QQMPolyRingElem};
    savelibname::String=joinpath(homedir(), "libs"),
    maple_path::String="maple",
    msolve_path::String="msolve",
    info_level::Int=0
)::Vector{Vector{Vector{QQFieldElem}}}
    procedure = """
        local vars, res:
        vars := [$(_get_str(gens(parent([eqs; pos; ineqs][1]))))]:
        res := RAG:-PointsPerComponents([$(_get_str(eqs))], [$(_get_str(pos))], [$(_get_str(ineqs))], {"verb"=$info_level}):
        return map(point -> subs(point, vars), res):
    """
    line = _call_procedure_with_maple(procedure; savelibname, maple_path, msolve_path)
    line = replace(line, r"(-?\d+)/(\d+)" => s"QQ(\1)/QQ(\2)")
    line = replace(line, r"(-?\d+)" => s"QQ(\1)")
    res = eval(Meta.parse(line))
    res
end

end # module RAG
