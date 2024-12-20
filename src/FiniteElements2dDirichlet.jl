module FiniteElements2dDirichlet

# Functions that the user can access
export init_LG_matrix, init_EQ_vector_and_m, init_Ke_matrix!, init_Fe_vector!, init_K_matrix, init_F_vector, init_mesh, solve_system, error_convergence, phi, f, d_xi_phi, x, d_xi_x

# Dependencies
using Plots
using GaussQuadrature
using SparseArrays
using LinearAlgebra


include("2d-stationary-dirichlet/phi_dphi.jl") # phi, d_xi_phi

include("2d-stationary-dirichlet/x_dx.jl") # x, d_xi_x

include("2d-stationary-dirichlet/LG.jl") # init_LG_matrix

include("2d-stationary-dirichlet/EQ_m.jl") # init_EQ_vector_and_m

include("2d-stationary-dirichlet/Ke_local.jl") # init_Ke_matrix!

include("2d-stationary-dirichlet/K_global.jl") # init_Fe_vector!

include("2d-stationary-dirichlet/Fe_local.jl") # init_Fe_vector!

include("2d-stationary-dirichlet/F_global.jl") # init_F_vector

include("2d-stationary-dirichlet/solve_system.jl") # solve_system

include("2d-stationary-dirichlet/error_convergence.jl") # error_convergence

include("2d-stationary-dirichlet/mesh.jl") # init_mesh

end
