"""
    init_F_vector(f, X_matrix, Y_matrix, m, EQ, LG)

Return a vector ``F`` with length `m` for solving the differential equation.

# Arguments
- `f::Function`: input function from the equation.
- `X_matrix::Matrix{Float64}`: x axis coordinates of the mesh generated by the init_mesh function.
- `Y_matrix::Matrix{Float64}`: y axis coordinates of the mesh generated by the init_mesh function.
- `m::Int`: m value generated by the EQ matrix generated by the init_EQ_vector_and_m function.
- `EQ::Matrix{Int}`: EQ matrix generated by the init_EQ_vector_and_m function.
- `LG::Matrix{Int}`: LG matrix generated by the init_LG_matrix function.

# Examples
```jldoctest
X, Y = init_mesh(4, 3);
LG = init_LG_matrix(4, 3);
EQ, m = init_EQ_vector_and_m(4, 3);
init_F_vector((x1, x2) -> 48, X, Y, m, EQ, LG);
# output
6-element Vector{Float64}:
 3.999999999999999
 3.999999999999999
 3.999999999999999
 3.999999999999999
 3.9999999999999996
 3.999999999999999
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