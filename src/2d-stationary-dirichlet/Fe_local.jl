"""
    init_Fe_vector!(f, Xs, Ys, Fe, P, W)

Mutates the given vector `Fe` (of length `4`) to compute the local `Fe` vector for finite element analysis.

# Arguments
- `f::Function`: input function from the equation.
- `Xs::Vector{Float64}`: vector of the element coordinates on the x axis.
- `Ys::Vector{Float64}`: vector of the element coordinates on the y axis.
- `Fe::Vector{Float64}`: vector that is mutated.
- `P::Vector{Float64}`: vector of gauss points for numerical integration.
- `W::Vector{Float64}`: vector of gauss weights for numerical integration.

# Examples
```jldoctest
using GaussQuadrature
Fe = zeros(4);
P, W = legendre(5);
init_Fe_vector!((x1, x2) -> 64, [0, 0.25, 0.25, 0], [0, 0, 0.25, 0.25], Fe, P, W);
Fe;
# output
4-element Vector{Float64}:
 0.9999999999999993
 0.9999999999999998
 1.0000000000000002
 0.9999999999999998
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