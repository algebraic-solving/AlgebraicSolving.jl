using AlgebraicSolving
using Documenter

push!(LOAD_PATH, "../src")

DocMeta.setdocmeta!(AlgebraicSolving, :DocTestSetup, :(using AlgebraicSolving); recursive=true)

# Build documentation.
# ====================

makedocs(
    # options
    modules = [AlgebraicSolving],
    doctest = true,
    clean = true,
    checkdocs = :none,
    sitename = "AlgebraicSolving.jl",
    format = Documenter.HTML(),
    pages = [
        "index.md",
        "types.md",
        "Algorithms" => ["groebner-bases.md",
                         "normal-forms.md",
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
