# Adds the path to the package to the load_path
push!(LOAD_PATH, "/home/joao_pereira/FiniteElements2dDirichlet.jl/src/")

using Documenter
using FiniteElements2dDirichlet

# Generates the documentation
makedocs(
    sitename = "FiniteElements2dDirichlet.jl",
)
