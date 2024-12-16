module FiniteElements2dDirichlet

# Functions that the user can access
export init_LG_matrix, init_EQ_vector_and_m, init_Ke_matrix!, init_Fe_vector!, init_K_matrix, init_F_vector, init_mesh, solve_system, error_convergence

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

"""
    init_LG_matrix(Nx, Ny)

Return the matrix that correlates the local to the global elements.

The size of the matrix is `4×Nx*Ny`

# Arguments
- `Nx::Integer`: the number of elements in the x axis.
- `Ny::Integer`: the number of elements in the y axis.

# Examples
```julia-repl
julia> init_LG_matrix(2, 4)
4×8 Matrix{Int64}:
1  2  4  5   7   8  10  11
2  3  5  6   8   9  11  12
5  6  8  9  11  12  14  15
4  5  7  8  10  11  13  14
```
"""
function init_LG_matrix(Nx, Ny)
  LG = fill(0, (4, Nx * Ny))
  j = 0

  for i in 1:Nx*Ny
    if (j % (Nx+1) == 0)
      j = j + 1
    end

    LG[1, i] = j
    LG[2, i] = j + 1
    LG[3, i] = j + Nx + 2
    LG[4, i] = j + Nx + 1

    j = j + 1
  end

  return LG
end


"""
    init_EQ_vector_and_m(Nx, Ny)

Return the vector that defines which terms are important and the dominium dimension.

The size of the vector is `Nx+1*Ny+1` and `m` is `(Nx-1)*(Ny-1)`.

# Arguments
- `Nx::Integer`: the number of elements in the x axis.
- `Ny::Integer`: the number of elements in the y axis.

# Examples
```julia-repl
julia> init_EQ_vector_and_m(2, 5)
([5, 5, 5, 5, 1, 5, 5, 2, 5, 5, 3, 5, 5, 4, 5, 5, 5, 5], 4)
```
"""
function init_EQ_vector_and_m(Nx, Ny)
  m = (Nx - 1) * (Ny - 1)
  EQ = zeros(Int, Ny+1, Nx+1)

  # Initialize the border elements
  for i in 1:Nx+1
    EQ[1,i] = m + 1
    EQ[Ny+1, i] = m + 1
  end
  for j in 1:Ny+1
    EQ[j, 1] = m+1
    EQ[j, Nx+1] = m+1
  end

  # initialize the within elements
  k = 1
  for i in 2:Ny
    for j in 2:Nx
      EQ[i,j] = k
      k = k + 1
    end
  end

  return cat(EQ'..., dims=1), m
end


"""
    init_Ke_matrix!(alpha, beta, Xs, Ys, Ke, P, W)

Mutates a given matrix `Ke` with dimensions `4×4` to become the local K matrix.

# Arguments
- `alpha::Float64`: constant α from the equation.
- `beta::Float64`: constant β from the equation.
- `Xs::Vector{Float64}`: vector of the element coordinates on the x axis.
- `Ys::Vector{Float64}`: vector of the element coordinates on the y axis.
- `Ke::Matrix{Float64}`: matrix that suffer the mutation.
- `P::Vector{Float64}`: vector of gauss points for numerical integration.
- `W::Vector{Float64}`: vector of gauss weights for numerical integration.

# Examples
```julia-repl
julia> init_Ke_matrix!(6, 0, [0, 0.25, 0.25, 0], [0, 0, 0.25, 0.25], Ke, P, W)
 4  -1  -2  -1
-1   4  -1  -2
-2  -1   4  -1
-1  -2  -1   4
```
"""
function init_Ke_matrix!(alpha, beta, Xs, Ys, Ke, P, W)
  fill!(Ke, 0)

  J_det = (xi1, xi2) -> abs(
    d_xi_x[1][1](xi1, xi2, Xs) * d_xi_x[2][2](xi1, xi2, Ys) -
    d_xi_x[1][2](xi1, xi2, Xs) * d_xi_x[2][1](xi1, xi2, Ys)
  )

  H = (xi1, xi2) -> [
     d_xi_x[2][2](xi1, xi2, Ys)   -d_xi_x[2][1](xi1, xi2, Ys);
    -d_xi_x[1][2](xi1, xi2, Xs)    d_xi_x[1][1](xi1, xi2, Xs)
  ]

  for i in 1:length(W)
    Wi, Pi = W[i], P[i]
    for j in 1:length(W)
      Wj, Pj = W[j], P[j]
      det_J = J_det(Pj, Pi)
      det_J_inv = 1 / det_J
      H_P = H(Pj, Pi)
      HTH = H_P' * H_P
      phi_vals = [phi[a](Pj, Pi) for a in 1:4]
      d_xi_phi_vals = [[d_xi_phi[1][a](Pj, Pi); d_xi_phi[2][a](Pj, Pi)] for a in 1:4]

      for a in 1:4
        grad_a = d_xi_phi_vals[a]
        aux = HTH * grad_a * det_J_inv * alpha

        for b in 1:4
          grad_b = d_xi_phi_vals[b]
          Ke[a, b] += Wi * Wj * ((grad_b' * aux) +  beta * phi_vals[b] * phi_vals[a] * det_J)
        end
      end
    end
  end
end


"""
    init_K_matrix(alpha, beta, X_matrix, Y_matrix, m, EQ, LG)

Return a matrix ``K`` with dimensions `m×m` for solving the differential equation.

# Arguments
- `alpha::Float64`: constant α from the equation.
- `beta::Float64`: constant β from the equation.
- `X_matrix::Matrix{Float64}`: x axis coordinates of the mesh generated by the init_mesh function.
- `Y_matrix::Matrix{Float64}`: y axis coordinates of the mesh generated by the init_mesh function.
- `m::Integer`: m value generated by the EQ matrix generated by the init_EQ_vector_and_m function.
- `EQ::Vector{Integer}`: EQ vector generated by the init_EQ_vector_and_m function.
- `LG::Matrix{Integer}`: LG matrix generated by the init_LG_matrix function.

# Examples
```julia-repl
julia> X, Y = init_mesh(2,2)
julia> EQ, m = init_EQ_vector_and_m(2,2)
julia> LG = init_LG_matrix(2,2)
julia> init_K_matrix(1, 1, X, Y, m, EQ, LG)
4×4 SparseArrays.SparseMatrixCSC{Float64, Int64} with 16 stored entries:
  2.71605   -0.320988  -0.320988  -0.330247
 -0.320988   2.71605   -0.330247  -0.320988
 -0.320988  -0.330247   2.71605   -0.320988
 -0.330247  -0.320988  -0.320988   2.71605
```
"""
function init_K_matrix(alpha, beta, X_matrix, Y_matrix, m, EQ, LG)
  ne = size(LG, 2) # assuming we are using a LG (4 x ne)
  K = spzeros(m+1, m+1)
  Ke = zeros(4, 4)
  P, W = legendre(5)

  for e in 1:ne
    init_Ke_matrix!(alpha, beta, X_matrix[LG[:,e]], Y_matrix[LG[:,e]], Ke, P, W)
    for b in 1:4
      j = EQ[LG[b, e]]
      for a in 1:4
        i = EQ[LG[a, e]]
        K[i,j] += Ke[a,b]
      end
    end
  end

  return K[1:m, 1:m]
end


"""
    init_Fe_vector!(f, Xs, Ys, Fe, P, W)

Mutates a given vector `Fe` with length `4` to become the local `Fe` vector.

# Arguments
- `f::Function`: input function from the equation.
- `Xs::Vector{Float64}`: vector of the element coordinates on the x axis.
- `Ys::Vector{Float64}`: vector of the element coordinates on the y axis.
- `Fe::Vector{Float64}`: vector that suffer the mutation.
- `P::Vector{Float64}`: vector of gauss points for numerical integration.
- `W::Vector{Float64}`: vector of gauss weights for numerical integration.

# Examples
```julia-repl
julia> Init_Fe_vector!((x1, x2) -> 64, [0, 0.25, 0.25, 0], [0, 0, 0.25, 0.25], Fe, P, W)
4-element Vector{Float64}:
1.0
1.0
1.0
1.0
```
"""
function init_Fe_vector!(f, Xs, Ys, Fe, P, W)
  fill!(Fe, 0)

  for i in 1:5
    Pi = P[i]
    for j in 1:5
      Pj = P[j]
      aux = f(x[1](Pi, Pj, Xs), x[2](Pi, Pj, Ys)) *
            abs(d_xi_x[1][1](Pi, Pj, Xs) *
                d_xi_x[2][2](Pi, Pj, Ys) -
                d_xi_x[1][2](Pi, Pj, Xs) *
                d_xi_x[2][1](Pi, Pj, Ys)) *
            W[i] * W[j]
      for a in 1:4
        Fe[a] += phi[a](Pi, Pj) * aux
      end
    end
  end
end


"""
    init_F_vector(f, X_matrix, Y_matrix, m, EQ, LG)

Return a vector ``F`` with length `m` for solving the differential equation.

# Arguments
- `f::Function`: input function from the equation.
- `X_matrix::Matrix{Float64}`: x axis coordinates of the mesh generated by the init_mesh function.
- `Y_matrix::Matrix{Float64}`: y axis coordinates of the mesh generated by the init_mesh function.
- `m::Integer`: m value generated by the EQ matrix generated by the init_EQ_vector_and_m function.
- `EQ::Matrix{Integer}`: EQ matrix generated by the init_EQ_vector_and_m function.
- `LG::Matrix{Integer}`: LG matrix generated by the init_LG_matrix function.

# Examples
```julia-repl
julia> X, Y = init_mesh(4, 3)
julia> LG = init_LG_matrix(4, 3)
julia> EQ, m = init_EQ_vector_and_m(4, 3)
julia> init_F_vector((x1, x2) -> 48, X, Y, m, EQ, LG)
6-element Vector{Float64}:
4.0
4.0
4.0
4.0
4.0
4.0
```
"""
function init_F_vector(f, X_matrix, Y_matrix, m, EQ, LG)
  ne = size(LG, 2) # assuming we are using a LG (4 x ne)
  F = zeros(m+1)
  Fe = zeros(4)
  P, W = legendre(5)

  for e in 1:ne
    init_Fe_vector!(f, X_matrix[LG[:,e]], Y_matrix[LG[:,e]], Fe, P, W)
    for a in 1:4
      F[EQ[LG[a,e]]] += Fe[a]
    end
  end

  return F[1: end-1]
end


"""
    solve_system(alpha, beta, f, Nx, Ny; EQLG=false, XY_matrix=false, noise=false)

Return the solution of the differential equation as a vector of the approximate solution
function evaluated in the inner knots of the mesh.

# Arguments
- `alpha::Float64`: constant α from the equation.
- `beta::Float64`: constant β from the equation.
- `f::Function`: input function from the equation.
- `Nx::Integer`: the number of elements in the x axis.
- `Ny::Integer`: the number of elements in the y axis.
- `EQLG::Bool`: indicates if the EQ vector and LG matrix should be returned.
- `XY_matrix::Bool`: indicates if the X and Y matrices from the mesh should be returned.
- `noise::Bool`: indicates if noise should be added to the mesh.

# Examples
```julia-repl
julia> solve_system(1, 1, (x, y) -> (2*pi^2 +1)*sin(pi*x)*sin(pi*y), 2, 3)
2-element Vector{Float64}:
 1.0040408191040564
 1.0040408191040564
```
"""
function solve_system(alpha, beta, f, Nx, Ny; EQLG=false, XY_matrix=false, noise=false)
  X_matrix, Y_matrix = init_mesh(Nx, Ny, ns=noise)
  EQ, m = init_EQ_vector_and_m(Nx, Ny)
  LG = init_LG_matrix(Nx, Ny)
  K = init_K_matrix(alpha, beta, X_matrix, Y_matrix, m, EQ, LG)
  F = init_F_vector(f, X_matrix, Y_matrix, m, EQ, LG)
  C = K \ F
  return EQLG ? XY_matrix ? (C, EQ, LG, X_matrix, Y_matrix) : (C, EQ, LG) : XY_matrix ? (C, X_matrix, Y_matrix) : C
end


# # Plot the approximation for a linear base
# function plot_approximation(alpha, beta, f, Nx1, Nx2; ns=false)
#
#   C, X, Y = solve_system(alpha, beta, f, Nx1, Nx2, XY_matrix=true, noise=ns)
#
#   # Initialize the axes and computes the temperatures (including the boundary condition) (works only in linear bases)
#   C_in_plane = vcat(zeros(1, Nx2+1), hcat(zeros(Nx1-1, 1), reshape(C, Nx1-1, Nx2-1), zeros(Nx1-1, 1)), zeros(1, Nx2+1))'
#
#   # Plots the approximation
#   # wireframe(X', Y', C_in_plane, linecolor=:black, lw=1, n=5, size=(500,500))
#   surface(X', Y', C_in_plane, color=:thermal, alpha=0.1, title="Approximation found for u(x1, x2)", xlabel="x1", ylabel="x2", zlabel="Temperatures", n=5)
#
#   savefig("approximation_found.png")
# end


"""
    error_convergence(lb, ub, alpha, beta, u, f; see_plot=false, ns=false)

return the errors of a given solution for discretizations with knots ranging from `2^lb` to `2^ub`.
The number of knots increases in powers of 2, starting from `2^lb` up to `2^ub`.

# Arguments
- `lb::Integer`: lower-bound limit for testing the error convergence.
- `ub::Integer`: upper-bound limit for testing the error convergence.
- `alpha::Float64`: constant α from the equation.
- `beta::Float64`: constant β from the equation.
- `u::Function`: solution of the equation.
- `f::Function`: input function from the equation.
- `see_plot::Bool`: indicates if the error convergence analysis should be plotted.
- `ns::Bool`: indicates if noise should be added to the mesh.

# Examples
```julia-repl
julia> error_convergence(2, 4, 1, 1, (x,y) -> sin(pi * x) * sin(pi * y), (x,y) -> (2*pi^2 + 1) * sin(pi * x) * sin(pi * y))
3-element Vector{Float64}:
 0.02946436125614044
 0.007348142187147988
 0.001836025960324437
```
"""
function error_convergence(lb, ub, alpha, beta, u, f; see_plot=false, ns=false)

  # Compute the error of a system given C
  function gauss_error(u, C, X_matrix, Y_matrix, EQ, LG)
    C_ext = [C; 0]  # Extend C with a zero
    erro = 0.0
    ne = size(LG, 2)  # Number of elements
    P, W = legendre(5)  # Gauss points and weights

    for e in 1:ne
      Xs, Ys = X_matrix[LG[:, e]], Y_matrix[LG[:, e]]
      element_error = 0.0

      for i in 1:5, j in 1:5
        xi1, xi2 = P[i], P[j]
        Wi, Wj = W[i], W[j]

        # Compute shape functions and their contributions
        x1 = x[1](xi1, xi2, Xs)
        x2 = x[2](xi1, xi2, Ys)
        aprox = sum(C_ext[EQ[LG[k, e]]] * phi[k](xi1, xi2) for k in 1:4)

        # Compute determinant of the Jacobian
        J_det = d_xi_x[1][1](xi1, xi2, Xs) * d_xi_x[2][2](xi1, xi2, Ys) -
                d_xi_x[1][2](xi1, xi2, Xs) * d_xi_x[2][1](xi1, xi2, Ys)

        # Integrand calculation
        element_error += Wi * Wj * (u(x1, x2) - aprox)^2 * J_det
      end

      erro += element_error
    end

    return sqrt(erro)
  end

  Nxs = [2^i for i in lb:ub]
  errors = zeros(ub - lb + 1)

  for i in lb:ub
    C, EQ, LG, X_matrix, Y_matrix = solve_system(alpha, beta, f, Nxs[i-lb+1], Nxs[i-lb+1]; EQLG=true, XY_matrix=true, noise=ns)
    errors[i-lb+1] = gauss_error(u, C, X_matrix , Y_matrix, EQ, LG)
  end

  if see_plot
    # Plot the errors in a log scale
    hs = 1 ./ Nxs
    plot(hs, errors, seriestype = :scatter, label = "Error convergence",
    xlabel = "h", ylabel = "error", size=(800, 800), xscale=:log10, yscale=:log10,
    markercolor = :blue)
    plot!(hs, errors, seriestype = :line, label = "", linewidth = 2, linecolor = :blue)
    plot!(hs, hs.^2, seriestype = :line, label = "h^2", linewidth = 2, linecolor = :red)

    # Save the graph in a png file
    savefig("error-convergence.png")
  end
  return errors
end


"""
    init_mesh(Nx, Ny; ns=false, plot=false)

Return the x and y axis mesh for solving the system.

# Arguments
- `Nx::Integer`: the number of elements in the x axis.
- `Ny::Integer`: the number of elements in the y axis.
- `ns::Bool`: indicates if noise should be added to the mesh.
- `plot::Bool`: indicates if the mesh should be plotted.

# Examples
```julia-repl
julia> init_mesh(2, 3)
([0.0 0.0 0.0 0.0; 0.5 0.5 0.5 0.5; 1.0 1.0 1.0 1.0], [0.0 0.3333333333333333 0.6666666666666666 1.0; 0.0 0.3333333333333333 0.6666666666666666 1.0; 0.0 0.3333333333333333 0.6666666666666666 1.0])
```
"""
function init_mesh(Nx, Ny; ns=false, plot=false)

  x1 = collect(0: 1/Nx :1)
  x2 = collect(0: 1/Ny :1)

  X = [x1[i] for i in 1:Nx+1, j in 1:Ny+1]
  Y = [x2[j] for i in 1:Nx+1, j in 1:Ny+1]

  if ns == true
    @assert size(X) == size(Y) "X e Y must have the same dimensions."
    if size(X, 1) > 2 && size(X, 2) > 2
      rl1, rl2 = 1 / 4*Nx, 1 / 4*Ny
      X[2:end-1, 2:end-1] .+= rl1 * (rand(Float64, size(X[2:end-1,2:end-1])) .- 0.5) * 2
      Y[2:end-1, 2:end-1] .+= rl2 * (rand(Float64, size(Y[2:end-1,2:end-1])) .- 0.5) * 2
    end
  end

  if plot == true
    Plots.plot(legend=false, aspect_ratio=:equal, xticks=0:0.25:1, yticks=0:0.25:1)
    Plots.scatter!(X, Y, markersize=4, color=:blue)

    LG = init_LG_matrix(Nx, Ny)
    ne = size(LG, 2)
    for e in 1:ne
      i1, i2, i3, i4 = LG[:, e]
      Plots.plot!([X[i1], X[i2], X[i3], X[i4], X[i1]],
                  [Y[i1], Y[i2], Y[i3], Y[i4], Y[i1]],
                  color=:black)
    end

    savefig("2d-mesh.png")
  end

  return X, Y
end

end
