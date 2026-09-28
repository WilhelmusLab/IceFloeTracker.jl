@testitem "LopezAcosta2019.Segment – simple case" tags = [:e2e, :smoke] begin
    import DataFrames: DataFrame, nrow
    import Images: RGB
    dataset = Watkins2026Dataset(; ref="v0.2")

    case = first(filter(c -> (c.case_number == 6 && c.satellite == "terra"), dataset))
    segments = LopezAcosta2019.Segment()(
        RGB.(modis_truecolor(case)), RGB.(modis_falsecolor(case)), modis_landmask(case)
    )
    expected_segment_count = validated_floe_properties(case) |> DataFrame |> nrow
    @test length(segments.segment_labels) ≈ expected_segment_count rtol = 0.7
end

@testitem "LopezAcosta2019.Segment – sample of cases" setup = [Segmentation] tags = [:e2e] begin
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
        LopezAcosta2019.Segment();
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
    @test mean_recall ≥ 0.8 broken = true
    @test mean_precision ≥ 0.8 broken = true
    @test mean_F_score ≥ 0.8 broken = true

    # Current performance should look at least as good as this:
    @test mean_recall ≥ 0.28
    @test mean_precision ≥ 0.3
    @test round(mean_F_score; digits=1) ≥ 0.38

    # return current performance
    @show mean_recall
    @show mean_precision
    @show mean_F_score
end

@testitem "LopezAcosta2019.Segment – detailed tests" setup = [Segmentation] tags = [:e2e] begin
    dataset = Watkins2026Dataset(; ref="v0.2")
    (; labeled_fraction, recall, precision, F_score) = run_and_validate_segmentation(
        first(filter(c -> (c.case_number == 6 && c.satellite == "terra"), dataset)),
        LopezAcosta2019.Segment();
        output_directory="./test_outputs/",
    )

    @test 0.10 ≈ labeled_fraction atol = 0.1
    @test 0.32 ≤ round(recall; digits=2)
    @test 0.92 ≤ round(precision; digits=2)
    @test 0.48 ≤ round(F_score; digits=2)

    (; labeled_fraction, recall, precision, F_score) = run_and_validate_segmentation(
        first(filter(c -> (c.case_number == 14 && c.satellite == "aqua"), dataset)),
        LopezAcosta2019.Segment();
        output_directory="./test_outputs/",
    )

    @test 0.05 ≈ labeled_fraction atol = 0.1 # lowered to 0.115
    @test 0.34 ≤ round(recall; digits=2)
    @test 0.85 ≤ round(precision; digits=2) # Note: Decreased precision, I suspect an issue with Seg. A.
    @test 0.49 ≤ round(F_score; digits=2)

    (; labeled_fraction, recall, precision, F_score) = run_and_validate_segmentation(
        first(filter(c -> (c.case_number == 61 && c.satellite == "aqua"), dataset)),
        LopezAcosta2019.Segment();
        output_directory="./test_outputs/",
    )
    
    @test 0.11 ≈ labeled_fraction atol = 0.1 # lowered to 0.21
    @test 0.30 ≤ round(recall; digits=2) # lowered to 0.54
    @test 0.96 ≤ round(precision; digits=2)
    @test 0.46 ≤ round(F_score; digits=2)

    (; labeled_fraction, recall, precision, F_score) = run_and_validate_segmentation(
        first(filter(c -> (c.case_number == 63 && c.satellite == "aqua"), dataset)),
        LopezAcosta2019.Segment();
        output_directory="./test_outputs/",
    )
    
    @test labeled_fraction ≈ 0.38 rtol = 0.1
    @test 0.60 ≤ round(recall; digits=2)
    @test 0.99 ≤ round(precision; digits=2)
    @test 0.75 ≤ round(F_score; digits=2)
end

@testitem "LopezAcosta2019.Segment – image types" setup = [Segmentation] tags = [:e2e] begin
    import Images: RGB, RGBA, n0f8, n6f10, n4f12, n2f14, n0f16, float32, float64
    dataset = Watkins2026Dataset(; ref="v0.2")
    case::Case = first(filter(c -> (c.case_number == 6 && c.satellite == "aqua"), dataset))
    algorithm = LopezAcosta2019.Segment()
    baseline = run_and_validate_segmentation(
        case, algorithm; output_directory="./test_outputs/"
    )

    @test results_invariant_for(RGB; baseline, algorithm, case)
    @test results_invariant_for(RGBA; baseline, algorithm, case)
    @test results_invariant_for(n0f8; baseline, algorithm, case)
    @test results_invariant_for(n6f10; baseline, algorithm, case) broken = true
    @test results_invariant_for(n4f12; baseline, algorithm, case) # broken = true
    @test results_invariant_for(n2f14; baseline, algorithm, case) # broken = true
    @test results_invariant_for(n0f16; baseline, algorithm, case)
    @test results_invariant_for(float32; baseline, algorithm, case)
    @test results_invariant_for(float64; baseline, algorithm, case)
    @test results_invariant_for(RGB, n0f8; baseline, algorithm, case)
    @test results_invariant_for(RGB, n6f10; baseline, algorithm, case) broken = true
    @test results_invariant_for(RGB, n4f12; baseline, algorithm, case) # broken = true
    @test results_invariant_for(RGB, n2f14; baseline, algorithm, case) # broken = true
    @test results_invariant_for(RGB, n0f16; baseline, algorithm, case)
    @test results_invariant_for(RGB, float32; baseline, algorithm, case)
    @test results_invariant_for(RGB, float64; baseline, algorithm, case)
end
