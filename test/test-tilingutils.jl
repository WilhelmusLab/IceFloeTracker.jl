using IceFloeTracker: get_optimal_tile_size
gots = get_optimal_tile_size
@testset "Tiling utils" begin

    @test gots(2, (10, 10)) == 2 # disregard tiles of 1 pixel
    @test gots(3, (15, 15)) == 3
    @test gots(4, (20, 20)) == 5 # prefer larger tile size

    # Test with edge case for minimum l0
    @test gots(3, (5, 5)) == 4
    @test gots(5, (5, 5)) == 5

    # Test with non-square dimensions
    @test gots(3, (10, 20)) == 2
    @test gots(3, (10, 7)) == 2

    # Test error handling for invalid l0
    @test_throws ErrorException gots(1, (10, 10))
    @test_throws ErrorException gots(7, (5, 5)) == 5
end
