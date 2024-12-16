using FiniteElements2dDirichlet
using LinearAlgebra
using Test
using GaussQuadrature

@testset "FiniteElements2dDirichlet.jl" begin

#---LG-EQ-m---
    # Testing if LG is being initialized correctly (this test is not attached to the Nx variables)
    @test init_LG_matrix(2, 4) == [1 2 4 5 7 8 10 11; 2 3 5 6 8 9 11 12; 5 6 8 9 11 12 14 15; 4 5 7 8 10 11 13 14]

    # Testing if EQ and m are being initializes correctly (this test is not attached to the Nx variables)
    @test init_EQ_vector_and_m(2, 5) == ([5, 5, 5, 5, 1, 5, 5, 2, 5, 5, 3, 5, 5, 4, 5, 5, 5, 5], 4)
#---LG-EQ-m---

#---Ke---
    # Initializing Ke
    Ke = zeros(4,4)
    P, W = legendre(5)

    # Testing the alpha component of the local K matrix initializer
    init_Ke_matrix!(6, 0, [0, 0.25, 0.25, 0], [0, 0, 0.25, 0.25], Ke, P, W)
    @test norm(Ke - [4 -1 -2 -1; -1 4 -1 -2; -2 -1 4 -1; -1 -2 -1 4]) < 0.001

    # Testing the beta component of the local K matrix initializer
    init_Ke_matrix!(0, 576.0, [0, 0.25, 0.25, 0.0], [0, 0, 0.25, 0.25], Ke, P, W)
    @test norm(Ke - [4 2 1 2; 2 4 2 1; 1 2 4 2; 2 1 2 4]) < 0.001

    # Testing both components of the local K matrix initializer
    init_Ke_matrix!(1, 1, [0, 2, 3, 1], [0, 0, 1, 1], Ke, P, W)
    @test norm(Ke - [0.72 0.11 0.06 -0.39; 0.11 1.72 -0.39 -0.94; 0.06 -0.39 0.72 0.11; -0.39 -0.94 0.11 1.72]) < 0.1
#---Ke---

#---K---
    X, Y = init_mesh(3,3)
    EQ, m = init_EQ_vector_and_m(3,3)
    LG = init_LG_matrix(3,3)
    @test norm(init_K_matrix(1, 1, X, Y, m, EQ, LG) - [2.71 -0.32 -0.32 -0.33; -0.32 2.71 -0.33 -0.32; -0.32 -0.33 2.72 -0.32; -0.33 -0.32 -0.32 2.72]) < 0.1
#---K---

#---Fe---
    # Initializing Fe
    Fe = zeros(4)
    P, W = legendre(5)

    # Testing the alpha component of the local F vector initializer
    init_Fe_vector!((x1, x2) -> 64, [0, 0.25, 0.25, 0], [0, 0, 0.25, 0.25], Fe, P, W)
    @test norm(Fe - [1; 1; 1; 1]) < 0.001

    # Testing the beta component of the local F vector initializer
    init_Fe_vector!((x1, x2) -> 36864*x1*x2, [0, 0.25, 0.25, 0.0], [0, 0, 0.25, 0.25], Fe, P, W)
    @test norm(Fe - [4; 8; 16; 8]) < 0.001

    # Testing both components of the local F vector initializer
    init_Fe_vector!((x1, x2) -> x1 + x2, [0, 2, 3, 1], [0, 0, 1, 1], Fe, P, W)
    @test norm(Fe - [2/3; 1; 4/3; 1]) < 0.001
#---Fe---

#---F---
    X, Y = init_mesh(4, 3)
    LG = init_LG_matrix(4, 3)
    EQ, m = init_EQ_vector_and_m(4, 3)
    @test norm(init_F_vector((x1, x2) -> 48, X, Y, m, EQ, LG) - [4; 4; 4; 4; 4; 4]) < 0.001

    X, Y = init_mesh(4, 4)
    LG = init_LG_matrix(4, 4)
    EQ, m = init_EQ_vector_and_m(4, 4)
    @test norm(init_F_vector((x1, x2) -> 36864*x1*x2, X, Y, m, EQ, LG) - [144; 288; 432; 288; 576; 864; 432; 864; 1296]) < 0.001

    # Testing with an irregular mesh
    X = [0 0 0 0; 0.25 0.266168 0.275391 0.25; 0.5 0.493792 0.521668 0.5; 0.75 0.747176 0.708237 0.75; 1 1 1 1]
    Y = [0 0.333333 0.666667 1; 0 0.352246 0.633228 1; 0 0.36139  0.693524 1; 0 0.326172 0.689905 1;  0 0.333333 0.666667 1]
    LG = init_LG_matrix(4, 3)
    EQ, m = init_EQ_vector_and_m(4, 3)
    @test norm(init_F_vector((x1, x2) -> x1 + x2, X, Y, m, EQ, LG) - [0.048; 0.069; 0.094; 0.077; 0.086; 0.116]) < 0.01
#---F---


#---init-mesh---
    X, Y = init_mesh(2, 3)
    @test norm(X - [0 0 0 0; 0.5 0.5 0.5 0.5; 1 1 1 1]) < 0.001
    @test norm(Y - [0 1/3 2/3 1.0; 0 1/3 2/3 1; 0 1/3 2/3 1]) < 0.001
#---init-mesh---

#---error-convergence---
    # Bounds for testing
    lb = 2
    ub = 4

    # Constants
    alpha = 1
    beta = 1

    # Functions
    u = (x1,x2) -> sin(pi * x1) * sin(pi * x2)
    ux1x1 = (x1,x2) -> -1 * pi^2 * sin(pi * x1) * sin(pi * x2)
    ux2x2 = (x1,x2) -> -1 * pi^2 * sin(pi * x1) * sin(pi * x2)
    f = (x1,x2) -> (-1 * alpha * ux1x1(x1,x2)) +
                   (-1 * alpha * ux2x2(x1,x2)) +
                   ( 1 * beta * u(x1,x2))

    #
    hs = 1 ./ [2^i for i in lb:ub]
    errors = error_convergence(lb, ub, alpha, beta, u, f)

    # Testing
    @test all(x -> abs(x - 2) < 0.1,
              map((i) -> log(errors[i+1] / errors[i]) /
                         log(hs[i+1] / hs[i]),
                  1:length(errors)-1))
#---error-convergence---

#---solve_system---
    # Mesh and Discretization
    Nx1 = 2
    Nx2 = 3
    X, Y = init_mesh(Nx1, Nx2)

    # Constants
    alpha = 1
    beta = 1

    # Functions
    u     = (x1,x2) -> sin(pi * x1) * sin(pi * x2)
    ux1x1 = (x1,x2) -> -1 * pi^2 * sin(pi * x1) * sin(pi * x2)
    ux2x2 = (x1,x2) -> -1 * pi^2 * sin(pi * x1) * sin(pi * x2)
    f     = (x1,x2) -> (-1 * alpha * ux1x1(x1,x2)) +
                       (-1 * alpha * ux2x2(x1,x2)) +
                       ( 1 * beta * u(x1,x2))

    # Testing
    @test norm(solve_system(alpha, beta, f, Nx1, Nx2) - vec(u.(X, Y)[2:end-1, 2:end-1])) < 1
#---solve_system---

end
