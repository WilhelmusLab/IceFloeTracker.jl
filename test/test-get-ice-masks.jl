
@testitem "get_ice_labels_mask tests" begin
    using IceFloeTracker:
        get_ice_labels_mask,
        get_tiles,
        kmeans_segmentation,
        get_nlabel,
        get_ice_masks,
        get_ice_peaks,
        apply_landmask

    include("config.jl")

    begin
        # Counts from case "111-greenland_sea-100km-20120623-terra-250m"
        edges = 0.0:0.015625:0.984375
        counts = [0.22849, 0.04548, 0.0317, 0.02294, 0.01893, 0.01466, 0.01293, 0.011,
                  0.00977, 0.00934, 0.00841, 0.00838, 0.00746, 0.00743, 0.00736, 0.00682,
                  0.00684, 0.00667, 0.00638, 0.00712, 0.00615, 0.0065, 0.00614, 0.00601,
                  0.00612, 0.00683, 0.00734, 0.00782, 0.00882, 0.0092, 0.00828, 0.00926,
                  0.00923, 0.00933, 0.00934, 0.00978, 0.0096, 0.01093, 0.01162, 0.01295,
                  0.0149, 0.01803, 0.02223, 0.0239, 0.02755, 0.03081, 0.03796, 0.04513,
                  0.054, 0.06292, 0.06747, 0.07155, 0.08288, 0.07175, 0.09395, 0.0995,
                  0.00321, 0.00075, 0.00022, 4.0e-5, 1.0e-5, 0.0, 0.0, 0.0]
        maximum_loc = get_ice_peaks(edges, counts)
        print(maximum_loc)
        @test maximum_loc == 0.859375
    end


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

