
using StatsBase: mean

@ntestset "preprocess_tiling" begin
    data_loader = Watkins2025GitHub(; ref="a451cd5e62a10309a9640fbbe6b32a236fcebc70")
    @ntestset "Detailed checks" begin
        (; labeled_fraction, recall, precision, F_score) = run_and_validate_segmentation(
            first(data_loader(c -> (c.case_number == 6 && c.satellite == "terra"))),
            LopezAcosta2019Tiling();
        )
        @test 0.426 ≈ labeled_fraction atol = 0.1
        @test 0.876 ≤ recall
        @test 0.595 ≤ precision
        @test 0.708 ≤ F_score

        (; labeled_fraction, recall, precision, F_score) = run_and_validate_segmentation(
            first(data_loader(c -> (c.case_number == 14 && c.satellite == "aqua"))),
            LopezAcosta2019Tiling();
        )
        @test 0.334 ≈ labeled_fraction atol = 0.1
        @test 0.846 ≤ recall
        @test 0.313 ≤ precision
        @test 0.457 ≤ F_score

        (; labeled_fraction, recall, precision, F_score) = run_and_validate_segmentation(
            first(data_loader(c -> (c.case_number == 61 && c.satellite == "aqua"))),
            LopezAcosta2019Tiling();
        )
        @test 0.271 ≈ labeled_fraction atol = 0.1
        @test 0.709 ≤ recall
        @test 0.686 ≤ precision
        @test 0.697 ≤ F_score

        (; labeled_fraction, recall, precision, F_score) = run_and_validate_segmentation(
            first(data_loader(c -> (c.case_number == 63 && c.satellite == "aqua"))),
            LopezAcosta2019Tiling();
        )
        @test 0.579 ≈ labeled_fraction atol = 0.1
        @test 0.901 ≤ recall
        @test 0.620 ≤ precision
        @test 0.734 ≤ F_score
    end

    @ntestset "Aggregate results" begin
        results = run_and_validate_segmentation(
            data_loader(case -> (case.case_number % 17 == 0)),
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
    @ntestset "Image types" begin
        case::ValidationDataCase = first(
            data_loader(c -> (c.case_number == 6 && c.satellite == "terra"))
        )
        algorithm = LopezAcosta2019Tiling()
        baseline = run_and_validate_segmentation(
            case, algorithm; output_directory="./test_outputs/"
        )
        function results_invariant_under_image_type_transformation(target_type::Function)
            intermediate_results_callback = save_results_callback(
                "./test_outputs/segmentation-LopezAcosta2019Tiling-$(target_type)-$(Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS"))";
            )

            segments = LopezAcosta2019Tiling()(
                target_type.(RGB.(case.modis_truecolor)),
                target_type.(RGB.(case.modis_falsecolor)),
                target_type.(RGB.(case.modis_landmask));
                intermediate_results_callback,
            )
            (; segment_count, labeled_fraction) = segmentation_summary(segments)
            (; recall, precision, F_score) = segmentation_comparison(
                case.validated_labeled_floes, segments
            )
            return all([
                ≈(segment_count, baseline.segment_count; rtol=0.01),
                ≈(labeled_fraction, baseline.labeled_fraction; rtol=0.01),
                ≈(recall, baseline.recall; rtol=0.01),
                ≈(precision, baseline.precision; rtol=0.01),
                ≈(F_score, baseline.F_score; rtol=0.01),
            ])
        end
        @test results_invariant_under_image_type_transformation(n0f8)
        @test results_invariant_under_image_type_transformation(n6f10) broken = true
        @test results_invariant_under_image_type_transformation(n4f12) broken = true
        @test results_invariant_under_image_type_transformation(n2f14) broken = true
        @test results_invariant_under_image_type_transformation(n0f16) broken = true
        @test results_invariant_under_image_type_transformation(float32) broken = true
        @test results_invariant_under_image_type_transformation(float64) broken = true
    end
end
