
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
        check_case(
            LopezAcosta2019Tiling(),
            data_loader,
            c -> (c.case_number == 6 && c.satellite == "terra");
            labeled_fraction_goal=0.426,
            recall_goal=0.876,
            precision_goal=0.595,
            F_score_goal=0.708,
        )
        check_case(
            LopezAcosta2019Tiling(),
            data_loader,
            c -> (c.case_number == 14 && c.satellite == "aqua");
            labeled_fraction_goal=0.334,
            recall_goal=0.846,
            precision_goal=0.313,
            F_score_goal=0.457,
        )
        check_case(
            LopezAcosta2019Tiling(),
            data_loader,
            c -> (c.case_number == 61 && c.satellite == "aqua");
            labeled_fraction_goal=0.271,
            recall_goal=0.709,
            precision_goal=0.686,
            F_score_goal=0.697,
        )
        check_case(
            LopezAcosta2019Tiling(),
            data_loader,
            c -> (c.case_number == 63 && c.satellite == "aqua");
            labeled_fraction_goal=0.579,
            recall_goal=0.901,
            precision_goal=0.620,
            F_score_goal=0.734,
        )
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
