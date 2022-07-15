using msolve_jll
using Nemo
using LinearAlgebra
using AlgebraicSolving
using Documenter

#= push!(LOAD_PATH, "../src") =#

include("../src/imports.jl")
include("../src/exports.jl")

DocMeta.setdocmeta!(AlgebraicSolving, :DocTestSetup, :(using AlgebraicSolving); recursive=true)

# Build documentation.
# ====================

makedocs(
    # options
    modules = [AlgebraicSolving],
    doctest = true,
    clean = true,
    strict = true,
    sitename = "AlgebraicSolving.jl",
    format = Documenter.HTML(),
    pages = [
        "index.md",
        "types.md",
        "Algorithms" => ["groebner-bases.md",
                         "solvers.md"],
        "Examples" => "katsura.md"
        ]
)

# Deploy built documentation from Travis.
# =======================================

deploydocs(
    # options
    repo = "github.com/algebraic-solving/AlgebraicSolving.jl.git",
    target = "build",
)
