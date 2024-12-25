"""
    init_EQ_vector_and_m(Nx, Ny)

Return the vector that defines which terms are important and the dominium dimension.

The size of the vector is `(Nx+1)*(Ny+1)` and `m` is `(Nx-1)*(Ny-1)`.

# Arguments
- `Nx::Integer`: the number of elements in the x axis.
- `Ny::Integer`: the number of elements in the y axis.

# Examples
```jldoctest
julia> init_EQ_vector_and_m(2, 5)
([5, 5, 5, 5, 1, 5, 5, 2, 5, 5, 3, 5, 5, 4, 5, 5, 5, 5], 4)
```
"""
function init_EQ_vector_and_m(Nx, Ny)
    @assert Nx > 0 "Nx must be a positive integer"
    @assert Ny > 0 "Ny must be a positive integer"
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