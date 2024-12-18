# Adiciona o caminho relativo ao diretório `src` no LOAD_PATH
push!(LOAD_PATH, "/home/joao_pereira/FiniteElements2dDirichlet.jl/src/")

using Documenter
using FiniteElements2dDirichlet

# Configura e gera a documentação
makedocs(
    sitename = "FiniteElements2dDirichlet.jl",
)
