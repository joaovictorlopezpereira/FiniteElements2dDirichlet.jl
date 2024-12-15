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

# Initializes the LG matrix
function init_LG_matrix(Nx1, Nx2)
  LG = fill(0, (4, Nx1 * Nx2))
  j = 0

  for i in 1:Nx1*Nx2
    if (j % (Nx1+1) == 0)
      j = j + 1
    end

    LG[1, i] = j
    LG[2, i] = j + 1
    LG[3, i] = j + Nx1 + 2
    LG[4, i] = j + Nx1 + 1

    j = j + 1
  end

  return LG
end

# Initializes the EQ Vector
function init_EQ_vector_and_m(Nx1, Nx2)
  m = (Nx1 - 1) * (Nx2 - 1)
  EQ = zeros(Int, Nx2+1, Nx1+1)

  # Initializes the border elements
  for i in 1:Nx1+1
    EQ[1,i] = m + 1
    EQ[Nx2+1, i] = m + 1
  end
  for j in 1:Nx2+1
    EQ[j, 1] = m+1
    EQ[j, Nx1+1] = m+1
  end

  # initializes the within elements
  k = 1
  for i in 2:Nx2
    for j in 2:Nx1
      EQ[i,j] = k
      k = k + 1
    end
  end

  return cat(EQ'..., dims=1), m
end

# Initializes the Ke matrix
function init_Ke_matrix!(alpha, beta, Xs, Ys, Ke, P, W)
  fill!(Ke, 0)

  # Determinante de Jacobiano
  J_det = (xi1, xi2) -> abs(
    d_xi_x[1][1](xi1, xi2, Xs) * d_xi_x[2][2](xi1, xi2, Ys) -
    d_xi_x[1][2](xi1, xi2, Xs) * d_xi_x[2][1](xi1, xi2, Ys)
  )

  # Matriz Jacobiana transposta vezes inversa
  H = (xi1, xi2) -> [
     d_xi_x[2][2](xi1, xi2, Ys)   -d_xi_x[2][1](xi1, xi2, Ys);
    -d_xi_x[1][2](xi1, xi2, Xs)    d_xi_x[1][1](xi1, xi2, Xs)
  ]

  # Loop pelos pontos de integração
  for i in 1:length(W)
    Wi, Pi = W[i], P[i]
    for j in 1:length(W)
      Wj, Pj = W[j], P[j]

      # Pré-computações locais
      det_J = J_det(Pj, Pi)
      det_J_inv = 1 / det_J
      H_P = H(Pj, Pi)
      HTH = H_P' * H_P

      # Pré-calcular `phi` e `d_xi_phi`
      phi_vals = [phi[a](Pj, Pi) for a in 1:4]
      d_xi_phi_vals = [[d_xi_phi[1][a](Pj, Pi); d_xi_phi[2][a](Pj, Pi)] for a in 1:4]

      for a in 1:4
        grad_a = d_xi_phi_vals[a]
        aux = HTH * grad_a * det_J_inv * alpha

        for b in 1:4
          grad_b = d_xi_phi_vals[b]
          Ke[a, b] += Wi * Wj * ((grad_b' * aux) +
                                 beta * phi_vals[b] * phi_vals[a] * det_J)
        end
      end
    end
  end
end

# Initializes the K matrix
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

# Initializes the Fe vector
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

# Initializes the F vector
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

# Solves the system using a regular mesh
function solve_system(alpha, beta, f, Nx1, Nx2; EQLG=false, XY_matrix=false, noise=false)
  X_matrix, Y_matrix = init_mesh(Nx1, Nx2, ns=noise)
  EQ, m = init_EQ_vector_and_m(Nx1, Nx2)
  LG = init_LG_matrix(Nx1, Nx2)
  K = init_K_matrix(alpha, beta, X_matrix, Y_matrix, m, EQ, LG)
  F = init_F_vector(f, X_matrix, Y_matrix, m, EQ, LG)
  C = K \ F
  return EQLG ? XY_matrix ? (C, EQ, LG, X_matrix, Y_matrix) : (C, EQ, LG) : XY_matrix ? (C, X_matrix, Y_matrix) : C
end

# # Plots the approximation for a linear base
# function plot_approximation(alpha, beta, f, Nx1, Nx2; ns=false)

#   C, X, Y = solve_system(alpha, beta, f, Nx1, Nx2, XY_matrix=true, noise=ns)

#   # Initializes the axes and computes the temperatures (including the boundary condition) (works only in linear bases)
#   C_in_plane = vcat(zeros(1, Nx2+1), hcat(zeros(Nx1-1, 1), reshape(C, Nx1-1, Nx2-1), zeros(Nx1-1, 1)), zeros(1, Nx2+1))'

#   # Plots the approximation
#   # wireframe(X', Y', C_in_plane, linecolor=:black, lw=1, n=5, size=(500,500))
#   surface(X', Y', C_in_plane, color=:thermal, alpha=0.1, title="Approximation found for u(x1, x2)", xlabel="x1", ylabel="x2", zlabel="Temperatures", n=5)

#   savefig("approximation_found.png")
# end

# Plots the error converge
function error_convergence(lb, ub, alpha, beta, u, f; see_plot=false, ns=false)

  # Computes the error of a system given C
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
    # Plots the errors in a log scale
    hs = 1 ./ Nxs
    plot(hs, errors, seriestype = :scatter, label = "Error convergence",
    xlabel = "h", ylabel = "error", size=(800, 800), xscale=:log10, yscale=:log10,
    markercolor = :blue)
    plot!(hs, errors, seriestype = :line, label = "", linewidth = 2, linecolor = :blue)
    plot!(hs, hs.^2, seriestype = :line, label = "h^2", linewidth = 2, linecolor = :red)

    # Saves the graph in a png file
    savefig("error-convergence.png")
  end
  return errors
end

# Initializes Xs and Ys as a regular mesh
function init_mesh(Nx1, Nx2; ns=false, plot=false)

  x1 = collect(0: 1/Nx1 :1)
  x2 = collect(0: 1/Nx2 :1)

  X = [x1[i] for i in 1:Nx1+1, j in 1:Nx2+1]
  Y = [x2[j] for i in 1:Nx1+1, j in 1:Nx2+1]

  if ns == true
    @assert size(X) == size(Y) "X e Y must have the same dimensions"
    if size(X, 1) > 2 && size(X, 2) > 2
      rl1, rl2 = 1 / 4*Nx1, 1 / 4*Nx2
      X[2:end-1, 2:end-1] .+= rl1 * (rand(Float64, size(X[2:end-1,2:end-1])) .- 0.5) * 2
      Y[2:end-1, 2:end-1] .+= rl2 * (rand(Float64, size(Y[2:end-1,2:end-1])) .- 0.5) * 2
    end
  end

  if plot == true
    Plots.plot(legend=false, aspect_ratio=:equal, xticks=0:0.25:1, yticks=0:0.25:1)
    Plots.scatter!(X, Y, markersize=4, color=:blue)

    LG = init_LG_matrix(Nx1, Nx2)
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
