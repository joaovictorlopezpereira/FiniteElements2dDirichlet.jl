module FiniteElements2dDirichlet

# Functions that the user can access
export init_LG_matrix, init_EQ_vector_and_m, init_Ke_matrix!, init_Fe_vector!, init_K_matrix, init_F_vector, init_mesh, solve_system, error_convergence, phi, f, d_xi_phi, x, d_xi_x

# Dependencies
using Plots
using GaussQuadrature
using SparseArrays
using LinearAlgebra


include("phi_dphi.jl") # phi, d_xi_phi

include("x_dx.jl") # x, d_xi_x

include("LG.jl") # init_LG_matrix

include("EQ_m.jl") # init_EQ_vector_and_m

include("Ke_local.jl") # init_Ke_matrix!

include("K_global.jl") # init_Fe_vector!

include("Fe_local.jl") # init_Fe_vector!

include("F_global.jl") # init_F_vector

include("solve_system.jl") # solve_system

include("error_convergence.jl") # error_convergence

include("mesh.jl") # init_mesh

end
