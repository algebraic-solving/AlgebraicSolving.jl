#
# parse arguments
#
length(ARGS) >= 1 || error("must provide path of msolve override directory as first argument")
msolveoverride = popfirst!(ARGS)

isdir(msolveoverride) || error("The given override path '$(msolveoverride)' is not a valid directory")
msolveoverride = abspath(msolveoverride)

#
#
#
@info "Install needed packages"
using Pkg
Pkg.develop(path=dirname(@__DIR__))
Pkg.add(["msolve_jll"])
Pkg.instantiate()

#
#
#
function add_jll_override(depot, pkgname, newdir)
    pkgid = Base.identify_package("$(pkgname)_jll")
    pkguuid = string(pkgid.uuid)
    mkpath(joinpath(depot, "artifacts"))
    open(joinpath(depot, "artifacts", "Overrides.toml"), "a") do f
        write(f, """
        [$(pkguuid)]
        $(pkgname) = "$(newdir)"
        """)
    end

    # we need to make sure that precompilation is run again with the override in place
    # (just running Pkg.precompile() does not seem to suffice)
    run(`touch $(Base.locate_package(pkgid))`)
end

tmpdepot = mktempdir(; cleanup=true)
@info "Created temporary depot at $(tmpdepot)"

# create override file for msolve_jll
add_jll_override(tmpdepot, "msolve", msolveoverride)

# create a fresh marker file in deps/src so the tree hash changes
#libmsolve_src_dir = joinpath(dirname(@__DIR__), "deps", "src")
#isdir(libmsolve_src_dir) || error("Could not find $(libmsolve_src_dir)")
#
#marker_prefix = ".recompile-libmsolve-julia-"
#for filename in readdir(libmsolve_src_dir)
#    if startswith(filename, marker_prefix)
#        rm(joinpath(libmsolve_src_dir, filename); force=true)
#    end
#end
#
#libmsolve_src_marker = joinpath(
#    libmsolve_src_dir,
#    marker_prefix * string(time_ns()) * ".txt",
#)
#write(libmsolve_src_marker, "force recompilation marker\n")

msolve_libdir = joinpath(msolveoverride, "lib")
dyld_fallback = let existing = get(ENV, "DYLD_FALLBACK_LIBRARY_PATH", "")
    isempty(existing) ? msolve_libdir : existing * ":" * msolve_libdir
end

# prepend our temporary depot to the depot list...
try
    withenv(
        "JULIA_DEPOT_PATH"=>tmpdepot*":"*join(DEPOT_PATH, ":"),
        "DYLD_FALLBACK_LIBRARY_PATH"=>dyld_fallback,
    ) do

        # ... and start Julia, by default with the same project environment
        run(`$(Base.julia_cmd()) --project=$(Base.active_project()) $(ARGS)`)
    end
finally
   # rm(libmsolve_src_marker; force=true)
end
