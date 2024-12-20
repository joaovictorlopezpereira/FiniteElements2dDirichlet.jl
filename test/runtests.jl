using FiniteElements2dDirichlet
using LinearAlgebra
using Test
using GaussQuadrature
using Documenter

DocMeta.setdocmeta!(FiniteElements2dDirichlet, :DocTestSetup, :(using FiniteElements2dDirichlet); recursive=true)

#---doctests---
doctest(FiniteElements2dDirichlet; manual=false)
#---doctests---

#---tests---
include("LG_tests.jl")
include("EQ_m_tests.jl")
include("Ke_tests.jl")
include("K_tests.jl")
include("Fe_tests..jl")
include("F_tests.jl")
include("mesh_tests.jl")
include("error_convergence_tests.jl")
include("solve_system_tests.jl")
#---tests---
