
@testitem "tiled ice labels" begin
    using IceFloeTracker:
        get_ice_peaks, get_tiles, get_ice_labels_mask, get_nlabel_relaxation
    using Images: build_histogram
    using Random
    include("config.jl")

    @testset "get_ice_peaks" begin
        Random.seed!(123)
        img = Gray.(rand(0:255, 10, 10) ./ 255)
        edges, counts = build_histogram(img, 64; minval=0, maxval=1)
        pk = get_ice_peaks(edges, counts)

        @test pk == 0.375
    end

    ref_img = load(falsecolor_test_image_file)
    tiles = get_tiles(ref_img; rblocks=8, cblocks=6)
    tile = tiles[1]
    thresholds = [10, 118, 120] ./ 255
    morph_residue = Gray.(readdlm("test_inputs/morph_residue_tile.csv", ',', Int) ./ 255)

    @testset "get_ice_labels" begin
        # regular use case applies landmask
        @test sum(get_ice_labels_mask(ref_img[tile...], thresholds)) == 6515
    end

    @testset "get_nlabel_relaxation" begin
        # regular use case applies landmask
        @test get_nlabel_relaxation(
            ref_img[tile...], morph_residue[tile...], 75/255, 10/255, 230/255
        ) == -1 # In the original, it defaulted to 1 if no ice was found.
    end
end