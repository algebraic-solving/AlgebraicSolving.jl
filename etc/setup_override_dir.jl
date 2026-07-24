# This Julia script sets up a directory with msolve compiled from a local path,
# then "installs" this msolve for use by the `run_with_override.jl` script

#
# parse arguments
#
length(ARGS) >= 1 || error("must provide path of msolve source directory as first argument")
length(ARGS) >= 2 || error("must provide path of destination directory as second argument")
msolve_prefix = popfirst!(ARGS)
prefix = popfirst!(ARGS)

if length(ARGS) > 0 && !startswith(ARGS[1], "--")
  build_dir = popfirst!(ARGS)
else
  build_dir = mktempdir(; cleanup = true)
end

run_configure = true
overwrite_allow = false
verbose = false
left_ARGS = String[]
while !isempty(ARGS)
   arg = popfirst!(ARGS)
   if arg == "--no-configure"
      global run_configure = false
   elseif arg == "--yes"
      global overwrite_allow = true
   elseif arg == "--verbose"
      global verbose = true
   else
      push!(left_ARGS, arg)
   end
end


# validate arguments
isdir(msolve_prefix) || error("The given msolve prefix '$(msolve_prefix)' is not a valid directory")
if ispath(prefix)
   if !overwrite_allow
      print("The given installation prefix '$(prefix)' already exists. Overwrite? [Y/n] ")
      overwrite_allow = lowercase(readline()) in ["y", "yes", ""]
   end
   if overwrite_allow
      rm(prefix; force=true, recursive=true)
   else
      error("Aborting")
   end
end

# convert into absolute paths
mkpath(prefix)
prefix = abspath(prefix)
msolve_prefix = abspath(msolve_prefix)
mkpath(build_dir)
build_dir = abspath(build_dir)

#
# Install needed packages
#
@info "Install needed packages"
using Pkg
using Artifacts
Pkg.add(["JLLPrefixes"])
Pkg.instantiate()

using JLLPrefixes

cd(build_dir)

if run_configure
   @info "Configuring msolve in $(build_dir) for $(prefix)"

   deps = ["GMP_jll", "FLINT_jll", "MPFR_jll"]
   artifact_paths = collect_artifact_paths(deps)
   deps_path = mktempdir(; cleanup=false)
   deploy_artifact_paths(deps_path, artifact_paths) # collect all (transitive) dependencies into one tree

   extraargs = [
        "CFLAGS=-I$(joinpath(deps_path, "include"))",
        "LDFLAGS=-L$(joinpath(deps_path, "lib"))",
   ]

   configure_cmd = `$(msolve_prefix)/configure
       --prefix=$(prefix)
       $(extraargs)
       $(left_ARGS)
       `

   verbose && @show configure_cmd

   # TODO: redirect the output of configure into a log file
   @show run(configure_cmd)
end


@info "Building msolve in $(build_dir)"
run(`make -j$(Sys.CPU_THREADS) $(verbose ? "V=1" : [])`)

@info "Installing msolve to $(prefix)"
run(`make install $(verbose ? "V=1" : [])`)
