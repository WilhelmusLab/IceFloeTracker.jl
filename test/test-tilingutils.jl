using IceFloeTracker:
    get_optimal_tile_size, get_tile_meta, bump_tile, get_tiles, to_uint8, get_brighten_mask, imbrighten
gots = get_optimal_tile_size

@testset "Tiling utils" begin
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

    @testset "bump_tile" begin
        extrarows, extracols = rand(1:100, 2)
        bumpby = (extrarows, extracols)
        @test bump_tile(tile, bumpby) == (1:(2 + extrarows), 3:(4 + extracols))
    end

    @testset "get_tiles" begin
        # unadjusted tiles
        _get_tiles(array, side_length) =
            collect(TileIterator(axes(array), (side_length, side_length)))

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

    @testset "get_brighten_mask" begin
        img = rand(0:255, 5, 5)
        bumped_img = img .+ 1
        mask = get_brighten_mask(img, bumped_img)
        @test all(mask .== 0)
    end

    @testset "imbrighten tests" begin
        img = [1 2; 3 4]
        brighten_mask = [1 0; 1 0]

        test_cases = [(1.25, [1 2; 4 4]), (0.1, [0 2; 0 4]), (0.9, img)]

        for (bright_factor, expected_result) in test_cases
            result = imbrighten(img, brighten_mask, bright_factor)
            @test result == expected_result
        end
    end
end