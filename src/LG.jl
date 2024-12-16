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