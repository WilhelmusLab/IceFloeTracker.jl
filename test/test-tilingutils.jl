using IceFloeTracker: get_optimal_tile_size, get_tile_meta, bump_tile, adjust_edge_tiles
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


    tile = (1:2, 3:4)

    @testset "get_tile_meta" begin
        @test get_tile_meta(tile) == [1, 2, 3, 4]
    end

    @testset "bump_tile" begin
        extrarows, extracols = 1, 2
        bumpby = (extrarows, extracols)
        @test bump_tile(tile, bumpby) == (1:2+extrarows, 3:4+extracols)
    end

    @testset "adjust_edge_tiles" begin
        array = rand(40, 20)
        l = 6
        tile_size = (l, l)
        bumpby = mod.(size(array), l)

        tiles = TileIterator(axes(array), tile_size) |> collect

        adjusted_tiles = adjust_edge_tiles(deepcopy(tiles), bumpby)

        # Test adjusted_tiles have one fewer column
        @test all(size(tiles) .- size(adjusted_tiles) .== (0, 1))

        # Test right edge tiles are bumped correctly
        _, m, _, n = get_tile_meta(adjusted_tiles[end]) - get_tile_meta(tiles[end-1, end-1])
        @test (m, n) == bumpby
    end
end
