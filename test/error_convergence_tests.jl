@testset "Error Convergence Test" begin
    # Initialize lb, ub
    lb = 2
    ub = 4

    # Initialize alpha, beta
    alpha = 1
    beta = 1

    # Initialize u, f
    u = (x1,x2) -> sin(pi * x1) * sin(pi * x2)
    ux1x1 = (x1,x2) -> -1 * pi^2 * sin(pi * x1) * sin(pi * x2)
    ux2x2 = (x1,x2) -> -1 * pi^2 * sin(pi * x1) * sin(pi * x2)
    f = (x1,x2) -> (-1 * alpha * ux1x1(x1,x2)) +
                  (-1 * alpha * ux2x2(x1,x2)) +
                  ( 1 * beta * u(x1,x2))

    # Initialize hs
    hs = 1 ./ [2^i for i in lb:ub]

    # Initialize errors
    errors = error_convergence(lb, ub, alpha, beta, u, f)

    # Test the error plot slope
    @test all(x -> abs(x - 2) < 0.1,
              map((i) -> log(errors[i+1] / errors[i]) /
                        log(hs[i+1] / hs[i]),
                  1:length(errors)-1))
end