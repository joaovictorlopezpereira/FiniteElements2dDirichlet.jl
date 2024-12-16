
@testset "K Initializer         " begin
    # Initialize X, Y, EQ, m, LG
    X, Y = init_mesh(3,3)
    EQ, m = init_EQ_vector_and_m(3,3)
    LG = init_LG_matrix(3,3)

    # Test K
    @test norm(init_K_matrix(1, 1, X, Y, m, EQ, LG) - [2.71 -0.32 -0.32 -0.33; -0.32 2.71 -0.33 -0.32; -0.32 -0.33 2.72 -0.32; -0.33 -0.32 -0.32 2.72]) < 0.1
end