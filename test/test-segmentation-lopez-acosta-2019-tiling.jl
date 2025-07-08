
using IceFloeTracker: LopezAcosta2019Tiling

@ntestset "$(@__FILE__)" begin
    @ntestset "smoke test" begin
        truecolor = load(
            "./test_inputs/pipeline/input_pipeline/20220914.aqua.reflectance.250m.tiff"
        )

        falsecolor = load(
            "./test_inputs/pipeline/input_pipeline/20220914.aqua.reflectance.250m.tiff"
        )

        landmask = load("./test_inputs/pipeline/input_pipeline/landmask.tiff")

        segments = LopezAcosta2019Tiling()(truecolor, falsecolor, landmask)
        @show segments
        save(
            "./test_outputs/segmentation-Lopez-Acosta-2019-Tiling-smoketest-mean-labels" *
            Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
            ".png",
            map(i -> segment_mean(segments, i), labels_map(segments)),
        )
        @test length(segment_labels(segments)) == 92
    end

    @ntestset "image types" begin
        truecolor = load(
            "./test_inputs/pipeline/input_pipeline/20220914.aqua.reflectance.250m.tiff"
        )
        falsecolor = load(
            "./test_inputs/pipeline/input_pipeline/20220914.aqua.reflectance.250m.tiff"
        )
        landmask = load("./test_inputs/pipeline/input_pipeline/landmask.tiff")
        region = (200:400, 500:900)
        for target_type in [n0f8, n6f10, n4f12, n2f14, n0f16, float32, float64]
            @info "Image type: $target_type"
            segments = LopezAcosta2019Tiling(; tile_rblocks=1, tile_cblocks=2)(
                target_type.(truecolor[region...]),
                target_type.(falsecolor[region...]),
                target_type.(landmask[region...]),
            )
            @show segments
            save(
                "./test_outputs/segmentation-LopezAcosta2019Tiling-mean-labels_$(target_type)_$(Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS")).png",
                map(i -> segment_mean(segments, i), labels_map(segments)),
            )
            @test length(segment_labels(segments)) == 36
        end
    end

    @ntestset "detailed test" begin
        region = (1016:3045, 1486:3714)
        data_dir = joinpath(@__DIR__, "test_inputs")
        true_color_image = load(
            joinpath(data_dir, "beaufort-chukchi-seas_truecolor.2020162.aqua.250m.tiff")
        )
        ref_image = load(
            joinpath(data_dir, "beaufort-chukchi-seas_falsecolor.2020162.aqua.250m.tiff")
        )
        landmask = float64.(load(joinpath(data_dir, "matlab_landmask.png"))) .> 0

        # Crop images to region of interest
        true_color_image, ref_image, landmask = [
            img[region...] for img in (true_color_image, ref_image, landmask)
        ]

        algorithm = LopezAcosta2019Tiling(;
            tile_rblocks=2,
            tile_cblocks=3,

            # Ice labels thresholds
            ice_labels_prelim_threshold=110.0,
            ice_labels_band_7_threshold=200.0,
            ice_labels_band_2_threshold=190.0,
            ice_labels_ratio_lower=0.0,
            ice_labels_ratio_upper=0.75,
            r_offset=0.0,

            # Adaptive histogram equalization parameters
            adapthisteq_white_threshold=25.5,
            adapthisteq_entropy_threshold=4,
            adapthisteq_white_fraction_threshold=0.4,

            # Gamma parameters,
            gamma=1.5,
            gamma_factor=1.3,
            gamma_threshold=220.0,

            # Unsharp mask parameters,
            unsharp_mask_radius=10,
            unsharp_mask_amount=2.0,
            unsharp_mask_factor=255.0,

            # Brighten parameters,
            brighten_factor=0.1,

            # Preliminary ice mask parameters,
            prelim_icemask_radius=10,
            prelim_icemask_amount=2,
            prelim_icemask_factor=0.5,

            # Main ice mask parameters,
            icemask_band_7_threshold=5,
            icemask_band_2_threshold=230,
            icemask_band_1_threshold=240,
            icemask_band_7_threshold_relaxed=10,
            icemask_band_1_threshold_relaxed=190,
            icemask_possible_ice_threshold=75,
            icemask_n_clusters=3,
        )
        segments = algorithm(ref_image, true_color_image, landmask)
        @info segments
        @test length(segment_labels(segments)) == 1383
        save(
            "./test_outputs/segmentation-Lopez-Acosta-2019-Tiling-detailed-mean-labels" *
            Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
            ".png",
            map(i -> segment_mean(segments, i), labels_map(segments)),
        )
    end
end
