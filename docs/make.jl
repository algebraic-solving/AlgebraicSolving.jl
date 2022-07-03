using AbstractAlgebra
using Nemo
using AlgebraicSolving
using Documenter

push!(LOAD_PATH, "../src/")

DocMeta.setdocmeta!(AlgebraicSolving, :DocTestSetup, :(using AlgebraicSolving, AbstractAlgebra, Nemo); recursive=true)

# Build documentation.
# ====================

makedocs(
    # options
    modules = [AlgebraicSolving],
    doctest = true,
    clean = true,
    sitename = "AlgebraicSolving.jl",
    format = Documenter.HTML(),
    pages = Any[
        "Introduction" => "index.md",
        "User Guide" => Any[
            "Algorithms" => ["algorithms/groebner-bases.md",
                    "algorithms/solvers.md"],
            "Examples" => "examples/katsura.md"
        ],
    ],
    strict = true
)

# Deploy built documentation from Travis.
# =======================================

deploydocs(
    # options
    repo = "github.com/algebraic-solving/AlgebraicSolving.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
    devbranch = "main"
)
