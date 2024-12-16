@testset "Solve System Test     " begin
    # Initialize Nx, Ny, X, Y
    Nx = 2
    Ny = 3
    X, Y = init_mesh(Nx, Ny)

    # Initialize alpha, beta
    alpha = 1
    beta = 1

    # Initialize u, f
    u   = (x,y) -> sin(pi * x) * sin(pi * y)
    uxx = (x,y) -> -1 * pi^2 * sin(pi * x) * sin(pi * y)
    uyy = (x,y) -> -1 * pi^2 * sin(pi * x) * sin(pi * y)
    f   = (x,y) -> (-1 * alpha * uxx(x,y)) +
                   (-1 * alpha * uyy(x,y)) +
                   ( 1 * beta * u(x,y))

    # Test if the found solution is close enough to the exact one
    @test norm(solve_system(alpha, beta, f, Nx, Ny) - vec(u.(X, Y)[2:end-1, 2:end-1])) < 1
end