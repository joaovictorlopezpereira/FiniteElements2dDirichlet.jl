
@testset "EQ and m Initializer  " begin
    # Test EQ and m
    @test init_EQ_vector_and_m(2, 5) == ([5, 5, 5, 5, 1, 5, 5, 2, 5, 5, 3, 5, 5, 4, 5, 5, 5, 5], 4)
end