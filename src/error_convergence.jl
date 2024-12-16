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