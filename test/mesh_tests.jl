@testset "Mesh Initializer      " begin
    # Initialize X, Y
    X, Y = init_mesh(2, 3)

    # Test X
    @test norm(X - [0 0 0 0; 0.5 0.5 0.5 0.5; 1 1 1 1]) < 0.001

    # Test Y
    @test norm(Y - [0 1/3 2/3 1.0; 0 1/3 2/3 1; 0 1/3 2/3 1]) < 0.001
end