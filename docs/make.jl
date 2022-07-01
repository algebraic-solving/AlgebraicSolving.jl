using Documenter
using AlgebraicSolving
using AbstractAlgebra

DocMeta.setdocmeta!(AlgebraicSolving, :DocTestSetup, :(using AlgebraicSolving); recursive=true)

# Build documentation.
# ====================

makedocs(
    # options
    modules = [AlgebraicSolving],
    doctest = true,
    clean = true,
    sitename = "AlgebraicSolving.jl",
    format = Documenter.HTML(
        #= canonical = "https://juliadata.github.io/DataFrames.jl/stable/",
         = assets = ["assets/favicon.ico"],
         = edit_link = "main" =#
    ),
    pages = Any[
        "Introduction" => "index.md",
        "User Guide" => Any[
            "Algorithms" => "groebner-bases.md"
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
