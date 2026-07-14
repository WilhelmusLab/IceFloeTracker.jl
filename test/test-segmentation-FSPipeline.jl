@testitem "FSPipeline.Segment – simple case" tags = [:e2e, :smoke] begin
    import DataFrames: DataFrame, nrow
    import Images: RGB
    dataset = Watkins2026Dataset(; ref="v0.2")

    case = first(filter(c -> (c.case_number == 6 && c.satellite == "terra"), dataset))
    segments = FSPipeline.Segment()(
        RGB.(modis_truecolor(case)), RGB.(modis_falsecolor(case)), modis_landmask(case)
    )
    expected_segment_count = validated_floe_properties(case) |> DataFrame |> nrow
    @test length(segments.segment_labels) ≈ expected_segment_count rtol = 0.7
end

@testitem "FSPipeline.Segment – sample of cases" setup = [Segmentation] tags = [:e2e] begin
    import StatsBase: mean
    dataset = Watkins2026Dataset(; ref="v0.2")
    results = run_and_validate_segmentation(
        filter(
            c -> (
                c.visible_floes == "yes" &&
                c.cloud_fraction_manual < 0.5 &&
                c.case_number % 3 == 0
            ),
            dataset,
        ),
        FSPipeline.Segment();
        output_directory="./test_outputs/",
    )
    @test all(results.success)

    # Aggregate performance measures
    mean_recall = round(mean(skipnanormissing(results.recall)); digits=2)
    mean_precision = round(mean(skipnanormissing(results.precision)); digits=2)
    mean_F_score = round(mean(skipnanormissing(results.F_score)); digits=2)

    # Good performance might look liks this:
    @test mean_recall ≥ 0.9 broken = true
    @test mean_precision ≥ 0.9 broken = true
    @test mean_F_score ≥ 0.9 broken = true

    # Better performance might look like this:
    @test mean_recall ≥ 0.8 broken = true # Note: Increase in recall without increase in precision may indicate higher false-positive rate. 
    @test mean_precision ≥ 0.8 broken = true
    @test mean_F_score ≥ 0.8 broken = true

    # Current performance should look at least as good as this:
    @test mean_recall ≥ 0.6
    @test mean_precision ≥ 0.3
    @test round(mean_F_score; digits=1) ≥ 0.38

    # return current performance
    @show mean_recall
    @show mean_precision
    @show mean_F_score
end

@testitem "FSPipeline.Segment – detailed tests" setup = [Segmentation] tags = [:e2e] begin
    dataset = Watkins2026Dataset(; ref="v0.2")
    (; labeled_fraction, recall, precision, F_score) = run_and_validate_segmentation(
        first(filter(c -> (c.case_number == 6 && c.satellite == "terra"), dataset)),
        LopezAcosta2019.Segment();
        output_directory="./test_outputs/",
    )
   
    @test 0.36 ≈ labeled_fraction atol = 0.1
    @test 0.68 ≤ round(recall; digits=2)
    @test 0.59 ≤ round(precision; digits=2)
    @test 0.63 ≤ round(F_score; digits=2)

    (; labeled_fraction, recall, precision, F_score) = run_and_validate_segmentation(
        first(filter(c -> (c.case_number == 14 && c.satellite == "aqua"), dataset)),
        FSPipeline.Segment();
        output_directory="./test_outputs/",
    )
    @test 0.16 ≈ labeled_fraction atol = 0.1
    @test 0.79 ≤ round(recall; digits=2)
    @test 0.78 ≤ round(precision; digits=2)
    @test 0.80 ≤ round(F_score; digits=2)

    (; labeled_fraction, recall, precision, F_score) = run_and_validate_segmentation(
        first(filter(c -> (c.case_number == 61 && c.satellite == "aqua"), dataset)),
        FSPipeline.Segment();
        output_directory="./test_outputs/",
    )

    @test 0.13 ≈ labeled_fraction atol = 0.1
    @test 0.36 ≤ round(recall; digits=2)
    @test 0.86 ≤ round(precision; digits=2)
    @test 0.53 ≤ round(F_score; digits=2)

    (; labeled_fraction, recall, precision, F_score) = run_and_validate_segmentation(
        first(filter(c -> (c.case_number == 63 && c.satellite == "aqua"), dataset)),
        FSPipeline.Segment();
        output_directory="./test_outputs/",
    )
    # Note: Validation dataset currently doesn't include the floes intersecting the edge.
    # Improving the segmentation lowered the scores here due to these floes.
    
    @test labeled_fraction ≈ 0.56 rtol = 0.1
    @test 0.90 ≤ round(recall; digits=2)
    @test 0.98 ≤ round(precision; digits=2)
    @test 0.93 ≤ round(F_score; digits=2)
end

@testitem "FSPipeline.Segment – image types" setup = [Segmentation] tags = [:e2e] begin
    import Images: RGB, RGBA, n0f8, n6f10, n4f12, n2f14, n0f16, float32, float64
    dataset = Watkins2026Dataset(; ref="v0.2")
    case::Case = first(filter(c -> (c.case_number == 6 && c.satellite == "aqua"), dataset))
    algorithm = FSPipeline.Segment()
    baseline = run_and_validate_segmentation(
        case, algorithm; output_directory="./test_outputs/"
    )

    @test results_invariant_for(RGB; baseline, algorithm, case)
    @test results_invariant_for(RGBA; baseline, algorithm, case)
    @test results_invariant_for(n0f8; baseline, algorithm, case)
    @test results_invariant_for(n6f10; baseline, algorithm, case)
    @test results_invariant_for(n4f12; baseline, algorithm, case)
    @test results_invariant_for(n2f14; baseline, algorithm, case)
    @test results_invariant_for(n0f16; baseline, algorithm, case)
    @test results_invariant_for(float32; baseline, algorithm, case)
    @test results_invariant_for(float64; baseline, algorithm, case)
    @test results_invariant_for(RGB, n0f8; baseline, algorithm, case)
    @test results_invariant_for(RGB, n6f10; baseline, algorithm, case)
    @test results_invariant_for(RGB, n4f12; baseline, algorithm, case)
    @test results_invariant_for(RGB, n2f14; baseline, algorithm, case)
    @test results_invariant_for(RGB, n0f16; baseline, algorithm, case)
    @test results_invariant_for(RGB, float32; baseline, algorithm, case)
    @test results_invariant_for(RGB, float64; baseline, algorithm, case)
end
