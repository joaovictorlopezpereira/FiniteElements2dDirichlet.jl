"""
    init_mesh(Nx, Ny; ns=false, plot=false)

Return the x and y axis mesh for solving the system.

# Arguments
- `Nx::Int`: the number of elements in the x axis.
- `Ny::Int`: the number of elements in the y axis.
- `ns::Bool`: indicates if noise should be added to the mesh.
- `plot::Bool`: indicates if the mesh should be plotted.

# Examples
```jldoctest
julia> init_mesh(2, 3)
([0.0 0.0 0.0 0.0; 0.5 0.5 0.5 0.5; 1.0 1.0 1.0 1.0], [0.0 0.3333333333333333 0.6666666666666666 1.0; 0.0 0.3333333333333333 0.6666666666666666 1.0; 0.0 0.3333333333333333 0.6666666666666666 1.0])```
"""
function init_mesh(Nx, Ny; ns=false, plot=false)
    @assert Nx > 0 "Nx must be a positive integer"
    @assert Ny > 0 "Ny must be a positive integer"

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