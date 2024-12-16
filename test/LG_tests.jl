
@testset "LG Initializer        " begin
    # Test LG
    @test init_LG_matrix(2, 4) == [1 2 4 5 7 8 10 11; 2 3 5 6 8 9 11 12; 5 6 8 9 11 12 14 15; 4 5 7 8 10 11 13 14]
end