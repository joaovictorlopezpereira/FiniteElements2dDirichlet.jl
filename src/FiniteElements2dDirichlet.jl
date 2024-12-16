module FiniteElements2dDirichlet

# Functions that the user can access
export init_LG_matrix, init_EQ_vector_and_m, init_Ke_matrix!, init_Fe_vector!, init_K_matrix, init_F_vector, init_mesh, solve_system, error_convergence, phi, f, d_xi_phi, x, d_xi_x

# Dependencies
using Plots
using GaussQuadrature
using SparseArrays
using LinearAlgebra


# Phi function-vector
phi = [
  (xi1, xi2) -> (1 - xi1) * (1 - xi2) * (1 / 4);
  (xi1, xi2) -> (1 + xi1) * (1 - xi2) * (1 / 4);
  (xi1, xi2) -> (1 + xi1) * (1 + xi2) * (1 / 4);
  (xi1, xi2) -> (1 - xi1) * (1 + xi2) * (1 / 4)
]

# Derivative of the phi function-vector
d_xi_phi = [
  [((xi1, xi2) -> (-1 / 4) * (1 - xi2)),
   ((xi1, xi2) -> ( 1 / 4) * (1 - xi2)),
   ((xi1, xi2) -> ( 1 / 4) * (1 + xi2)),
   ((xi1, xi2) -> (-1 / 4) * (1 + xi2))],
  [((xi1, xi2) -> (-1 / 4) * (1 - xi1)),
   ((xi1, xi2) -> (-1 / 4) * (1 + xi1)),
   ((xi1, xi2) -> ( 1 / 4) * (1 + xi1)),
   ((xi1, xi2) -> ( 1 / 4) * (1 - xi1))]
]

# X function-vector
x = [
  (xi1, xi2, Xs) -> sum(Xs[k] * phi[k](xi1, xi2) for k in 1:4),
  (xi1, xi2, Ys) -> sum(Ys[k] * phi[k](xi1, xi2) for k in 1:4)
]

# Derivative of the x function-vector
d_xi_x = [
  [(xi1, xi2, coords) -> sum(coords[k] * d_xi_phi[1][k](xi1, xi2) for k in 1:4),
   (xi1, xi2, coords) -> sum(coords[k] * d_xi_phi[2][k](xi1, xi2) for k in 1:4)],
  [(xi1, xi2, coords) -> sum(coords[k] * d_xi_phi[1][k](xi1, xi2) for k in 1:4),
   (xi1, xi2, coords) -> sum(coords[k] * d_xi_phi[2][k](xi1, xi2) for k in 1:4)]
]

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
