
@testitem "FSPipeline.Segment – simple case" tags = [:e2e, :smoke] begin
    import DataFrames: DataFrame, nrow
    dataset = Watkins2026Dataset(; ref="v0.1")
    case = first(filter(c -> (c.case_number == 6 && c.satellite == "terra"), dataset))
    algo = FSPipeline.Segment()
    segments = FSPipeline.Segment()(
        modis_truecolor(case), modis_falsecolor(case), modis_landmask(case)
    )
    expected_segment_count = validated_floe_properties(case) |> DataFrame |> nrow
    @test length(segments.segment_labels) ≈ expected_segment_count rtol = 1.0
end
@testitem "FSPipeline.Segment - detailed" setup = [Segmentation] tags = [:e2e] begin
    dataset = Watkins2026Dataset(; ref="v0.1")
    (; labeled_fraction, recall, precision, F_score) = run_and_validate_segmentation(
        first(filter(c -> (c.case_number == 6 && c.satellite == "terra"), dataset)),
        FSPipeline.Segment();
    )
    @test 0.43 ≈ labeled_fraction atol = 0.1
    @test 0.87 ≤ round(recall; digits=2)
    @test 0.56 ≤ round(precision; digits=2)
    @test 0.69 ≤ round(F_score; digits=2)

    (; labeled_fraction, recall, precision, F_score) = run_and_validate_segmentation(
        first(filter(c -> (c.case_number == 14 && c.satellite == "aqua"), dataset)),
        FSPipeline.Segment();
    )
    @test 0.33 ≈ labeled_fraction atol = 0.1
    @test 0.66 ≤ round(recall; digits=2)
    @test 0.27 ≤ round(precision; digits=2)
    @test 0.38 ≤ round(F_score; digits=2)

    (; labeled_fraction, recall, precision, F_score) = run_and_validate_segmentation(
        first(filter(c -> (c.case_number == 61 && c.satellite == "aqua"), dataset)),
        FSPipeline.Segment();
    )
    @test 0.27 ≈ labeled_fraction atol = 0.1
    @test 0.71 ≤ round(recall; digits=2)
    @test 0.67 ≤ round(precision; digits=2)
    @test 0.70 ≤ round(F_score; digits=2)

    (; labeled_fraction, recall, precision, F_score) = run_and_validate_segmentation(
        first(filter(c -> (c.case_number == 63 && c.satellite == "aqua"), dataset)),
        FSPipeline.Segment();
    )
    @test 0.58 ≈ labeled_fraction atol = 0.1
    @test 0.83 ≤ round(recall; digits=2)
    @test 0.62 ≤ round(precision; digits=2)
    @test 0.73 ≤ round(F_score; digits=2)
end
# Segmentation fault with case 102 with Julia 1.12, passed with Julia 1.11
@testitem "FSPipeline.Segment - sample of cases" setup = [Segmentation] tags = [
    :e2e
] begin
    using StatsBase: mean
    dataset = Watkins2026Dataset(; ref="v0.1")

    #### Performance categories
    # 1. Clearly visible floes. The cases we want to see success with.
    # 2. No visible floes, either due to full cloud cover or resolution (e.g., slush, small floes, filaments)
    # 3. Partially visible floes -- ambiguous cases.

    results = run_and_validate_segmentation(
        filter(case -> (case.visible_floes == "yes" && case.cloud_fraction_manual <= 0.5 && case.case_number % 2 == 0), dataset),
        FSPipeline.Segment();
        output_directory="./test_outputs/",
    )
    @test all(results.success)

    results[:, ["case_name", "success", "recall", "precision", "F_score", "labeled_fraction", "segment_count"]] |> save("test_outputs/FSPipeline_sampled.csv")

    # Aggregate performance measures
    mean_recall = round(mean(skipnanormissing(results.recall)), digits=2)
    mean_precision = round(mean(skipnanormissing(results.precision)), digits=2)
    mean_F_score = round(mean(skipnanormissing(results.F_score)), digits=2)

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
    @test mean_precision ≥ 0.3
    @test round(mean_F_score; digits=1) ≥ 0.2

    # return current performance
    @show mean_recall
    @show mean_precision
    @show mean_F_score
end

@testitem "FSPipeline.Segment - image types" tags = [:e2e] setup = [Segmentation] begin
    using Images: RGB, RGBA, n0f8, n6f10, n4f12, n2f14, n0f16, float32, float64
    dataset = Watkins2026Dataset(; ref="v0.1")

    case = first(filter(c -> (c.case_number == 6 && c.satellite == "terra"), dataset))
    algorithm = FSPipeline.Segment()
    baseline = run_and_validate_segmentation(
        case, algorithm; output_directory="./test_outputs/"
    )

    @test results_invariant_for(RGB; baseline, algorithm, case)
    @test results_invariant_for(RGB, n0f8; baseline, algorithm, case)
    @test results_invariant_for(RGB, n6f10; baseline, algorithm, case)
    @test results_invariant_for(RGB, n4f12; baseline, algorithm, case)
    @test results_invariant_for(RGB, n2f14; baseline, algorithm, case)
    @test results_invariant_for(RGB, n0f16; baseline, algorithm, case)
    @test results_invariant_for(RGB, float32; baseline, algorithm, case)
    @test results_invariant_for(RGB, float64; baseline, algorithm, case)
    @test results_invariant_for(RGBA; baseline, algorithm, case)
    @test results_invariant_for(RGBA, n0f8; baseline, algorithm, case)
    @test results_invariant_for(RGBA, n6f10; baseline, algorithm, case)
    @test results_invariant_for(RGBA, n4f12; baseline, algorithm, case)
    @test results_invariant_for(RGBA, n2f14; baseline, algorithm, case)
    @test results_invariant_for(RGBA, n0f16; baseline, algorithm, case)
    @test results_invariant_for(RGBA, float32; baseline, algorithm, case)
    @test results_invariant_for(RGBA, float64; baseline, algorithm, case)
end
