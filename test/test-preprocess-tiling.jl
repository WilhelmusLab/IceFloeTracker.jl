
using IceFloeTracker:
    adapthisteq_params,
    adjust_gamma_params,
    brighten_factor,
    cloud_mask_thresholds,
    ice_masks_params,
    prelim_icemask_params,
    preprocess_tiling,
    structuring_elements,
    unsharp_mask_params,
    get_tiles,
    binarize_segments
using StatsBase: mean

skipnanormissing(arr::AbstractArray) = filter(x -> !ismissing(x) && !isnan(x), arr)

include("segmentation_utils.jl")

@testset "preprocess_tiling" begin
    result_images_to_save = [
        :ref_image,
        :true_color_image,
        :ref_img_cloudmasked,
        :prelim_icemask,
        :binarized_tiling,
        :segment_mask,
        :L0mask,
        :icemask,
        :final,
        :segment_mean_truecolor,
        :segment_mean_falsecolor,
    ]
    region = (1016:3045, 1486:3714)
    data_dir = joinpath(@__DIR__, "test_inputs")
    true_color_image = load(
        joinpath(data_dir, "beaufort-chukchi-seas_truecolor.2020162.aqua.250m.tiff")
    )
    ref_image = load(
        joinpath(data_dir, "beaufort-chukchi-seas_falsecolor.2020162.aqua.250m.tiff")
    )
    landmask_image = load(joinpath(data_dir, "matlab_landmask.png"))

    # Crop images to region of interest
    true_color_image, ref_image, landmask_image = [
        img[region...] for img in (true_color_image, ref_image, landmask_image)
    ]

    tile_settings = (; rblocks=2, cblocks=3)

    segments = LopezAcosta2019Tiling(;
        tile_settings,
        cloud_mask_thresholds,
        adapthisteq_params,
        adjust_gamma_params,
        structuring_elements,
        unsharp_mask_params,
        ice_masks_params,
        prelim_icemask_params,
        brighten_factor,
    )(
        true_color_image,
        ref_image,
        landmask_image;
        intermediate_results_callback=save_results_callback(
            "./test_outputs/segmentation-LopezAcosta2019Tiling-functor-$(Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS"))";
            names=result_images_to_save,
        ),
    )
    binary_floe_mask = binarize_segments(segments)
    @test abs(sum(binary_floe_mask) - 1461116) / 1461116 < 0.1

    # Same test, a different way
    (; labeled_fraction) = segmentation_summary(segments)
    @test labeled_fraction ≈ 0.3015 atol = 0.03

    @ntestset "Validated data" begin
        data_loader = Watkins2025GitHub(; ref="a451cd5e62a10309a9640fbbe6b32a236fcebc70")
        results = run_segmentation_over_multiple_cases(
            data_loader,
            case -> (case.case_number % 17 == 0),
            LopezAcosta2019Tiling();
            output_directory="./test_outputs/",
            result_images_to_save,
        )
        @test all(results.success)

        # Aggregate performance measures
        mean_recall = mean(skipnanormissing(results.recall))
        mean_precision = mean(skipnanormissing(results.precision))
        mean_F_score = mean(skipnanormissing(results.F_score))

        # Good performance might look liks this:
        @test mean_recall ≥ 0.9 broken = true
        @test mean_precision ≥ 0.9 broken = true
        @test mean_F_score ≥ 0.9 broken = true

        # Better performance might look like this:
        @test mean_recall ≥ 0.8 broken = true
        @test mean_precision ≥ 0.8 broken = true
        @test mean_F_score ≥ 0.8 broken = true

        # Current performance should look at least as good as this:
        @test mean_recall ≥ 0.5
        @test mean_precision ≥ 0.05
        @test mean_F_score ≥ 0.1

        # return current performance
        @info mean_recall
        @info mean_precision
        @info mean_F_score
    end
end
