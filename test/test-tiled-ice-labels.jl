
@testitem "tiled ice labels" setup = [Paths] begin
    using IceFloeTracker:
        get_image_peaks, get_tiles, get_ice_labels_mask, get_nlabel_relaxation
    using Random

    @testset "get_image_peaks" begin
        Random.seed!(123)
        img = rand(0:255, 10, 10)
        l, h = get_image_peaks(img)
        @test sum(l[1:5]) == 324
        @test sum(h[1:5]) == 11
    end

    ref_img = load(falsecolor_test_image_file)
    tiles = get_tiles(ref_img; rblocks=8, cblocks=6)
    tile = tiles[1]
    factor = 255
    thresholds = [10, 118, 120]
    morph_residue = readdlm("test_inputs/morph_residue_tile.csv", ',', Int)

    @testset "get_ice_labels" begin
        # regular use case applies landmask
        @test sum(get_ice_labels_mask(ref_img[tile...], thresholds, 255)) == 6515
    end

    @testset "get_nlabel_relaxation" begin
        # regular use case applies landmask
        @test get_nlabel_relaxation(
            ref_img[tile...], morph_residue[tile...], factor, 75, 10, 230
        ) == 1
    end
end