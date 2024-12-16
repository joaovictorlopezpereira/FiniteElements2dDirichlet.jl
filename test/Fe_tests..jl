
@testset "Fe Initializer        " begin
    # Initialize Fe, P, W
    Fe = zeros(4)
    P, W = legendre(5)

    # Test the alpha component of the local F vector initializer
    init_Fe_vector!((x1, x2) -> 64, [0, 0.25, 0.25, 0], [0, 0, 0.25, 0.25], Fe, P, W)
    @test norm(Fe - [1; 1; 1; 1]) < 0.001

    # Test the beta component of the local F vector initializer
    init_Fe_vector!((x1, x2) -> 36864*x1*x2, [0, 0.25, 0.25, 0.0], [0, 0, 0.25, 0.25], Fe, P, W)
    @test norm(Fe - [4; 8; 16; 8]) < 0.001

    # Test both components of the local F vector initializer
    init_Fe_vector!((x1, x2) -> x1 + x2, [0, 2, 3, 1], [0, 0, 1, 1], Fe, P, W)
    @test norm(Fe - [2/3; 1; 4/3; 1]) < 0.001
end