"""
    solve_system(alpha, beta, f, Nx, Ny; EQLG=false, XY_matrix=false, noise=false)

Return the solution of the differential equation as a vector of the approximate solution
function evaluated in the inner knots of the mesh.

# Arguments
- `alpha::Float64`: constant α from the equation.
- `beta::Float64`: constant β from the equation.
- `f::Function`: input function from the equation.
- `Nx::Int`: the number of elements in the x axis.
- `Ny::Int`: the number of elements in the y axis.
- `EQLG::Bool`: indicates if the EQ vector and LG matrix should be returned.
- `XY_matrix::Bool`: indicates if the X and Y matrices from the mesh should be returned.
- `noise::Bool`: indicates if noise should be added to the mesh.

# Examples
```jldoctest
julia> solve_system(1, 1, (x, y) -> (2*pi^2 +1)*sin(pi*x)*sin(pi*y), 2, 3)
2-element Vector{Float64}:
 1.0040408191040564
 1.0040408191040564
```
"""
function solve_system(alpha, beta, f, Nx, Ny; EQLG=false, XY_matrix=false, noise=false)
    @assert Nx > 0 "Nx must be a positive integer"
    @assert Ny > 0 "Ny must be a positive integer"

    X_matrix, Y_matrix = init_mesh(Nx, Ny, ns=noise)
    EQ, m = init_EQ_vector_and_m(Nx, Ny)
    LG = init_LG_matrix(Nx, Ny)
    K = init_K_matrix(alpha, beta, X_matrix, Y_matrix, m, EQ, LG)
    F = init_F_vector(f, X_matrix, Y_matrix, m, EQ, LG)
    C = K \ F
    return EQLG ? XY_matrix ? (C, EQ, LG, X_matrix, Y_matrix) : (C, EQ, LG) : XY_matrix ? (C, X_matrix, Y_matrix) : C
end