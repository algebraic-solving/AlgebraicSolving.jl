# How to use the override scripts

## Using AlgebraicSolving.jl with a different version of msolve than what `msolve_jll` provides

This can be useful for various reasons e.g.,

- you need to test msolve.jl with a newer msolve version, perhaps even its master branch,
- you need to test with a msolve with custom compiler settings (e.g. to enable debugging),
- and so on.

For this to work, follow these instructions:

1. Obtain a copy of the msolve sources, probably from a clone of the msolve git repository.
   Let's say this is in directory `MSOLVEROOT`.

2. Execute `./autogen.sh` in the `MSOLVEROOT` directory.

3. Build msolve by executing the `etc/setup_override_dir.jl` script in the `override` environment.
   Arguments:
    - first argument: the `MSOLVEROOT`
    - second argument: the directory where the override environment shall be installed, e.g. `/tmp/msolve_jll_override`.
    - third argument (optional): a temp build directory to make use of incremental builds,
      e.g. `/tmp/msolve_jll_override_build`. If not given, a temporary directory will be created
      and deleted after the build.
    - `--no-configure` (optional): if given, the script will not execute `./configure` in the `MSOLVEROOT`,
      but assume that this has already been done. This can be useful if you want to save time by doing incremental builds,
      or if you need to pass special arguments to `./configure` that the script does not know about.
    - `--yes` (optional): if given, the script will not ask for confirmation before deleting the override directory (if it already exists).
    - `--verbose` (optional): if given, the script will print more verbose output.

   To give a concrete example you could invoke

        julia --proj=override etc/setup_override_dir.jl $MSOLVEROOT /tmp/msolve_jll_override /tmp/msolve_jll_override_build

   for the first time, and then for subsequent builds (after making changes to the msolve sources) you could invoke

        julia --proj=override etc/setup_override_dir.jl $MSOLVEROOT /tmp/msolve_jll_override /tmp/msolve_jll_override_build --no-configure --yes

4. Use the `etc/run_with_override.jl` script with the exact same Julia executable
   and the override environment we just prepared.

        julia --proj=override etc/run_with_override.jl /tmp/msolve_jll_override

5. This opens a Julia session with the override in effect. You can now e.g. load AlgebraicSolving.jl
   via `using AlgebraicSolving`, or install other packages (such as Oscar) and test with them.
