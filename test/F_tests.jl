@testset "F Initializer         " begin
    # Initialize X, Y, LG, EQ, m
    X, Y = init_mesh(4, 3)
    LG = init_LG_matrix(4, 3)
    EQ, m = init_EQ_vector_and_m(4, 3)

    # Test F
    @test norm(init_F_vector((x1, x2) -> 48, X, Y, m, EQ, LG) - [4; 4; 4; 4; 4; 4]) < 0.001

    # Initialize X, Y, LG, EQ, m
    X, Y = init_mesh(4, 4)
    LG = init_LG_matrix(4, 4)
    EQ, m = init_EQ_vector_and_m(4, 4)

    # Test F
    @test norm(init_F_vector((x1, x2) -> 36864*x1*x2, X, Y, m, EQ, LG) - [144; 288; 432; 288; 576; 864; 432; 864; 1296]) < 0.001

    # Initialize X, Y, LG, EQ, m, with X and Y as irregular meshes
    X = [0 0 0 0; 0.25 0.266168 0.275391 0.25; 0.5 0.493792 0.521668 0.5; 0.75 0.747176 0.708237 0.75; 1 1 1 1]
    Y = [0 0.333333 0.666667 1; 0 0.352246 0.633228 1; 0 0.36139  0.693524 1; 0 0.326172 0.689905 1;  0 0.333333 0.666667 1]
    LG = init_LG_matrix(4, 3)
    EQ, m = init_EQ_vector_and_m(4, 3)

    # Test F
    @test norm(init_F_vector((x1, x2) -> x1 + x2, X, Y, m, EQ, LG) - [0.048; 0.069; 0.094; 0.077; 0.086; 0.116]) < 0.01
end