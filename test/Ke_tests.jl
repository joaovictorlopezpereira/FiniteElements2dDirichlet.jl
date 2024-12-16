
@testset "Ke Initializer        " begin
    # Initialize Ke
    Ke = zeros(4,4)
    P, W = legendre(5)

    # Test the alpha component of the local K matrix initializer
    init_Ke_matrix!(6, 0, [0, 0.25, 0.25, 0], [0, 0, 0.25, 0.25], Ke, P, W)
    @test norm(Ke - [4 -1 -2 -1; -1 4 -1 -2; -2 -1 4 -1; -1 -2 -1 4]) < 0.001

    # Test the beta component of the local K matrix initializer
    init_Ke_matrix!(0, 576.0, [0, 0.25, 0.25, 0.0], [0, 0, 0.25, 0.25], Ke, P, W)
    @test norm(Ke - [4 2 1 2; 2 4 2 1; 1 2 4 2; 2 1 2 4]) < 0.001

    # Test both components of the local K matrix initializer
    init_Ke_matrix!(1, 1, [0, 2, 3, 1], [0, 0, 1, 1], Ke, P, W)
    @test norm(Ke - [0.72 0.11 0.06 -0.39; 0.11 1.72 -0.39 -0.94; 0.06 -0.39 0.72 0.11; -0.39 -0.94 0.11 1.72]) < 0.1
end