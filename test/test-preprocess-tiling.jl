
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
    get_tiles

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
    @ntestset "established example" begin
        region = (1016:3045, 1486:3714)
        data_dir = joinpath(@__DIR__, "test_inputs")
        true_color_image = load(
            joinpath(data_dir, "beaufort-chukchi-seas_truecolor.2020162.aqua.250m.tiff")
        )
        ref_image = load(
            joinpath(data_dir, "beaufort-chukchi-seas_falsecolor.2020162.aqua.250m.tiff")
        )
        landmask_image = load(joinpath(data_dir, "matlab_landmask.png"))
        datestamp = Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS")

        # Crop images to region of interest
        true_color_image, ref_image, landmask_image = [
            img[region...] for img in (true_color_image, ref_image, landmask_image)
        ]

        tiles = get_tiles(true_color_image; rblocks=2, cblocks=3)
        @ntestset "function version" begin
            landmask = float64.(landmask_image) .> 0
            binary_floe_mask = preprocess_tiling(
                ref_image,
                true_color_image,
                (dilated=landmask,),
                tiles,
                cloud_mask_thresholds,
                adapthisteq_params,
                adjust_gamma_params,
                structuring_elements,
                unsharp_mask_params,
                ice_masks_params,
                prelim_icemask_params,
                brighten_factor;
                intermediate_results_callback=save_results_callback(
                    "./test_outputs/segmentation-LopezAcosta2019Tiling-function-$(datestamp)";
                    names=result_images_to_save,
                ),
            )

            # dmw: replace with test of mismatch against a preprocessed image
            @test abs(sum(binary_floe_mask) - 1461116) / 1461116 < 0.1

            labels = label_components(binary_floe_mask)
            segments = SegmentedImage(true_color_image, labels)
            (; labeled_fraction) = segmentation_summary(segments)
            @test labeled_fraction ≈ 0.3015 atol = 0.03
        end
        @ntestset "functor version" begin
            segments = LopezAcosta2019Tiling(;
                tile_settings=(; rblocks=2, cblocks=3),
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
                    "./test_outputs/segmentation-LopezAcosta2019Tiling-functor-$(datestamp)";
                    names=result_images_to_save,
                ),
            )
            binary_floe_mask = labels_map(segments) .> 0
            @test abs(sum(binary_floe_mask) - 1461116) / 1461116 < 0.1

            (; labeled_fraction) = segmentation_summary(segments)
            @test labeled_fraction ≈ 0.3015 atol = 0.03
        end
    end

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
    end
end
