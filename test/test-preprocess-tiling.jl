
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
    data_loader = Watkins2025GitHub(; ref="a451cd5e62a10309a9640fbbe6b32a236fcebc70")
    @ntestset "Detailed checks" begin
        @ntestset "Watkins 2025, case 14, aqua" begin
            case = first(
                data_loader(;
                    case_filter=c -> (c.case_number == 14 && c.satellite == "aqua")
                ),
            )

            validated_segments = case.validated_labeled_floes
            algorithm = LopezAcosta2019Tiling()
            measured_segments = algorithm(
                case.modis_truecolor,
                case.modis_falsecolor,
                case.modis_landmask;
                intermediate_results_callback=save_results_callback(
                    "./test_outputs/", case, algorithm
                ),
            )
            (; labeled_fraction) = segmentation_summary(measured_segments)
            (; precision, recall, F_score) = segmentation_comparison(
                validated_segments, measured_segments
            )
            @show labeled_fraction
            @show recall
            @show precision
            @show F_score

            @test labeled_fraction ≈ 0.3346 atol = 0.03
            @test 0.846 ≤ recall
            @test 0.313 ≤ precision
            @test 0.457 ≤ F_score
        end
    end
    @ntestset "Aggregate results" begin
        results = run_segmentation_over_multiple_cases(
            data_loader,
            case -> (case.case_number % 17 == 0),
            LopezAcosta2019Tiling();
            output_directory="./test_outputs/",
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
        @show mean_recall
        @show mean_precision
        @show mean_F_score
    end
end
