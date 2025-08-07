
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
    @ntestset "Validated data" begin
        data_loader = Watkins2025GitHub(; ref="a451cd5e62a10309a9640fbbe6b32a236fcebc70")
        @ntestset "visible floes, no clouds, no artifacts" begin
            dataset = data_loader(;
                case_filter=c -> (
                    c.visible_floes == "yes" &&
                    c.cloud_category_manual == "none" &&
                    c.artifacts == "no"
                ),
            )
            @info dataset.metadata

            results = run_segmentation_over_multiple_cases(
                dataset.data, LopezAcosta2019Tiling(); output_directory="./test_outputs/"
            )
            @info results

            # Run tests on aggregate results
            # First be sure we have the right number of results
            @test nrow(results) == nrow(dataset.metadata)

            # Now check that all cases run through without crashing
            successes = subset(results, :success => ByRow(==(true)))
            @test nrow(results) == nrow(successes)
        end

        @ntestset "visible floes, thin clouds, no artifacts" begin
            dataset = data_loader(;
                case_filter=c -> (
                    c.visible_floes == "yes" &&
                    c.cloud_category_manual == "thin" &&
                    c.artifacts == "no" &&
                    c.case_number % 5 == 0
                ),
            )
            @info dataset.metadata
            results = run_segmentation_over_multiple_cases(
                dataset.data, LopezAcosta2019Tiling(); output_directory="./test_outputs/"
            )
            @info results

            # Run tests on aggregate results
            # First be sure we have the right number of results
            @test nrow(results) == nrow(dataset.metadata)

            # Now check that all cases run through without crashing
            successes = subset(results, :success => ByRow(==(true)))
            @test nrow(results) == nrow(successes)
        end
        @ntestset "random sample" begin
            dataset = data_loader(; case_filter=c -> (c.case_number % 17 == 0))
            @info dataset.metadata
            results = run_segmentation_over_multiple_cases(
                dataset.data, LopezAcosta2019Tiling(); output_directory="./test_outputs/"
            )
            @info results

            # Run tests on aggregate results
            # First be sure we have the right number of results
            @test nrow(results) == nrow(dataset.metadata)

            # Now check that all cases run through without crashing
            successes = subset(results, :success => ByRow(==(true)))
            @test nrow(results) == nrow(successes)
        end
    end

    @ntestset "Full size" begin
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

        landmask = (dilated=landmask,)
        tiles = get_tiles(true_color_image; rblocks=2, cblocks=3)

        foo = preprocess_tiling(
            ref_image,
            true_color_image,
            landmask,
            tiles,
            cloud_mask_thresholds,
            adapthisteq_params,
            adjust_gamma_params,
            structuring_elements,
            unsharp_mask_params,
            ice_masks_params,
            prelim_icemask_params,
            brighten_factor,
        )

        # dmw: replace with test of mismatch against a preprocessed image
        @test abs(sum(foo) - 1461116) / 1461116 < 0.1
    end
end
