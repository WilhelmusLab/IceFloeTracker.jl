
@testitem "get_ice_labels_mask tests" begin
    using IceFloeTracker:
        get_ice_labels_mask,
        get_tiles,
        kmeans_segmentation,
        get_nlabel,
        get_ice_masks,
        apply_landmask

    include("config.jl")

    begin
        region = (1016:3045, 1486:3714)
        data_dir = joinpath(@__DIR__, "test_inputs")
        ref_image = load(
            joinpath(data_dir, "beaufort-chukchi-seas_falsecolor.2020162.aqua.250m.tiff")
        )
        landmask = float64.(load(joinpath(data_dir, "matlab_landmask.png"))) .> 0
        ref_image, landmask = [img[region...] for img in (ref_image, landmask)]
        morph_residue = readdlm(joinpath(data_dir, "ice_masks/morph_residue.csv"), ',', Int)
    end

    tiles = get_tiles(ref_image; rblocks=2, cblocks=3)
    ref_image_landmasked = apply_landmask(ref_image, .!landmask)

    begin
        tile = tiles[1]
        band_7_threshold = 5
        band_2_threshold = 230
        band_1_threshold = 240
        thresholds = (band_7_threshold, band_2_threshold, band_1_threshold)
        foo = get_ice_labels_mask(ref_image[tile...], thresholds)
        @test sum(foo) == 0

        morph_residue_seglabels = kmeans_segmentation(Gray.(morph_residue[tile...] / 255))
        @test get_nlabel(ref_image_landmasked[tile...], morph_residue_seglabels, 255) == 3
    end

    begin # first relaxation
        first_relax_thresholds = (10, band_2_threshold, 190)
        bar = get_ice_labels_mask(ref_image[tile...], first_relax_thresholds)
        @test sum(bar) == 8
    end

    begin
        tile = tiles[2]
        foo = get_ice_labels_mask(ref_image[tile...], thresholds)
        @test sum(foo) == 32
    end

    begin
        morph_residue_seglabels = kmeans_segmentation(Gray.(morph_residue[tile...] / 255))
        @test get_nlabel(ref_image_landmasked[tile...], morph_residue_seglabels, 255) == 3
    end

    begin
        tile = tiles[3]
        foo = get_ice_labels_mask(ref_image[tile...], thresholds)
        @test sum(foo) == 1
    end

    begin
        tile = tiles[4]
        foo = get_ice_labels_mask(ref_image[tile...], thresholds)
        @test sum(foo) == 29
    end

    begin
        tile = tiles[5]
        foo = get_ice_labels_mask(ref_image[tile...], thresholds)
        @test sum(foo) == 19
    end

    begin
        tile = tiles[6]
        foo = get_ice_labels_mask(ref_image[tile...], thresholds)
        @test sum(foo) == 62
    end

    begin
        morph_residue_seglabels = kmeans_segmentation(
            Gray.(morph_residue[tile...] / 255); k=3
        )
        @test get_nlabel(ref_image_landmasked[tile...], morph_residue_seglabels, 255) == 1
    end

    ice_mask, binarized_tiling = get_ice_masks(ref_image, morph_residue, landmask, tiles)
    @test sum(ice_mask) == 2669451
    @test sum(binarized_tiling) == 2873080
end