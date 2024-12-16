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
julia> Ke = zeros(4,4)
julia> P, W = legendre(5)
julia> init_Ke_matrix!(6, 0, [0, 0.25, 0.25, 0], [0, 0, 0.25, 0.25], Ke, P, W)
julia> Ke
4×4 Matrix{Float64}:
  4.0  -1.0  -2.0  -1.0
 -1.0   4.0  -1.0  -2.0
 -2.0  -1.0   4.0  -1.0
 -1.0  -2.0  -1.0   4.0
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