# Adiciona o caminho relativo ao diretório `src` no LOAD_PATH
push!(LOAD_PATH, "C:/Users/João/Desktop/FiniteElements2dDirichlet.jl-master")

using Documenter
using FiniteElements2dDirichlet

# Configura e gera a documentação
makedocs(
    sitename = "FiniteElements2dDirichlet.jl",
    remotes=nothing
)
