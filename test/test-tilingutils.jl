
@testitem "Tiling utils" begin
    using TiledIteration: TileIterator

    gots = get_optimal_tile_size

    @testset "get_optimal_tile_size" begin
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

    tile = (1:2, 3:4)

    @testset "get_tile_meta" begin
        @test get_tile_meta(tile) == [1, 2, 3, 4]
    end

    @testset "get_tiles" begin
        # unadjusted tiles
        _get_tiles(array, side_length) =
            TileIterator(axes(array), (side_length, side_length)) |> collect

        array = rand(40, 20)

        # Total coverage, no adjustment needed
        side_length = 5
        tiles = _get_tiles(array, side_length)
        @test tiles == get_tiles(array, side_length)

        # Leftovers greater than side_length ÷ 2, no adjustment needed
        side_length = 7 # bumpby = (5, 6)
        tiles = _get_tiles(array, side_length)
        adjusted_tiles = get_tiles(array, side_length)
        @test tiles == get_tiles(array, side_length)

        # adjust right edge tiles
        side_length = 6 # bumpby = (4, 2) => crop 1 column and bump by 2
        tiles = _get_tiles(array, side_length)
        adjusted_tiles = get_tiles(array, side_length)
        @test all(size(tiles) .- size(adjusted_tiles) .== (0, 1))

        # general case with both edges adjusted
        expected_tiles = [
            (1:5, 1:5) (1:5, 6:10) (1:5, 11:16)
            (6:12, 1:5) (6:12, 6:10) (6:12, 11:16)
        ]

        array = rand(12, 16)
        side_length = 5
        newtiles = get_tiles(array, side_length)
        @test expected_tiles == newtiles

        lowerleft_tile = newtiles[end, 1]
        _, b, _, _ = get_tile_meta(lowerleft_tile)

        lowerright_tile = newtiles[end, end]
        _, _, _, d = get_tile_meta(lowerright_tile)
        @test (b, d) == size(array)
    end
end

@testitem "get_tiles" begin
    @test get_tiles(rand(1, 1), 1) == [(1:1, 1:1);;]
    @test get_tiles(rand(1, 2), 1) == [(1:1, 1:1) (1:1, 2:2);]
    @test get_tiles(rand(1, 3), 1) == [(1:1, 1:1) (1:1, 2:2) (1:1, 3:3);]
    @test get_tiles(rand(2, 2), 1) == [(1:1, 1:1) (1:1, 2:2); (2:2, 1:1) (2:2, 2:2)]

    # These cases skip the last row of tiles.
    @test get_tiles(rand(4, 3), 2) == [(1:2, 1:3); (3:4, 1:3);;]
    @test get_tiles(rand(6, 4), 3) == [(1:3, 1:4); (4:6, 1:4);;]
    @test get_tiles(rand(10, 6), 5) == [(1:5, 1:6); (6:10, 1:6);;]
    @test get_tiles(rand(12, 16), 5) == [
        (1:5, 1:5) (1:5, 6:10) (1:5, 11:16)
        (6:12, 1:5) (6:12, 6:10) (6:12, 11:16)
    ]
    @test get_tiles(rand(20, 11), 10) == [(1:10, 1:11); (11:20, 1:11);;]
    @test get_tiles(rand(6000, 3556), (2000, 3556)) ==
        [(1:2000, 1:3556); (2001:4000, 1:3556); (4001:6000, 1:3556);;]
end

@testitem "MergeLastTileIfSmallerThanHalf" begin
    using TiledIteration: TileIterator
    using IceFloeTracker.ImageUtils: MergeLastTileIfSmallerThanHalf

    @test TileIterator((1:1,), MergeLastTileIfSmallerThanHalf((1,))) == [(1:1,)]
    @test TileIterator((1:2,), MergeLastTileIfSmallerThanHalf((1,))) == [(1:1,), (2:2,)]
    @test TileIterator((1:3,), MergeLastTileIfSmallerThanHalf((1,))) ==
        [(1:1,), (2:2,), (3:3,)]

    @test TileIterator((1:3,), MergeLastTileIfSmallerThanHalf((2,))) == [(1:2,), (3:3,)]
    @test TileIterator((1:4,), MergeLastTileIfSmallerThanHalf((2,))) == [(1:2,), (3:4,)]
    @test TileIterator((1:7,), MergeLastTileIfSmallerThanHalf((2,))) ==
        [(1:2,), (3:4,), (5:6,), (7:7,)]

    @test TileIterator((1:4,), MergeLastTileIfSmallerThanHalf((3,))) == [(1:4,)]
    @test TileIterator((1:5,), MergeLastTileIfSmallerThanHalf((3,))) == [(1:3,), (4:5,)]
    @test TileIterator((1:6,), MergeLastTileIfSmallerThanHalf((3,))) == [(1:3,), (4:6,)]
    @test TileIterator((1:7,), MergeLastTileIfSmallerThanHalf((3,))) == [(1:3,), (4:7,)]
    @test TileIterator((1:8,), MergeLastTileIfSmallerThanHalf((3,))) ==
        [(1:3,), (4:6,), (7:8,)]

    @test TileIterator((1:8,), MergeLastTileIfSmallerThanHalf((4,))) == [(1:4,), (5:8,)]
    @test TileIterator((1:9,), MergeLastTileIfSmallerThanHalf((4,))) == [(1:4,), (5:9,)]
    @test TileIterator((1:10,), MergeLastTileIfSmallerThanHalf((4,))) ==
        [(1:4,), (5:8,), (9:10,)]

    @test TileIterator((1:10,), MergeLastTileIfSmallerThanHalf((5,))) == [(1:5,), (6:10,)]
    @test TileIterator((1:12,), MergeLastTileIfSmallerThanHalf((5,))) == [(1:5,), (6:12,)]
    @test TileIterator((1:16,), MergeLastTileIfSmallerThanHalf((5,))) ==
        [(1:5,), (6:10,), (11:16,)]

    @test TileIterator((1:11,), MergeLastTileIfSmallerThanHalf((10,))) == [(1:11,)]

    @test TileIterator((1:6000,), MergeLastTileIfSmallerThanHalf((1000,))) ==
        [(1:1000,), (1001:2000,), (2001:3000,), (3001:4000,), (4001:5000,), (5001:6000,)]
    @test TileIterator((1:3556,), MergeLastTileIfSmallerThanHalf((1000,))) ==
        [(1:1000,), (1001:2000,), (2001:3000,), (3001:3556,)]
end