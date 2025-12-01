
@testitem "LopezAcosta2019Tiling.Segment – simple case" tags = [:e2e, :smoke] begin
    import DataFrames: DataFrame, nrow
    dataset = Watkins2026Dataset(; ref="v0.1")
    case = first(filter(c -> (c.case_number == 6 && c.satellite == "terra"), dataset))
    segments = LopezAcosta2019Tiling.Segment()(
        modis_truecolor(case), modis_falsecolor(case), modis_landmask(case)
    )
    expected_segment_count = validated_floe_properties(case) |> DataFrame |> nrow
    @test length(segments.segment_labels) ≈ expected_segment_count rtol = 1.0
end
@testitem "LopezAcosta2019Tiling.Segment - detailed" setup = [Segmentation] tags = [:e2e] begin
    dataset = Watkins2026Dataset(; ref="v0.1")
    (; labeled_fraction, recall, precision, F_score) = run_and_validate_segmentation(
        first(filter(c -> (c.case_number == 6 && c.satellite == "terra"), dataset)),
        LopezAcosta2019Tiling.Segment();
    )
    @test 0.43 ≈ labeled_fraction atol = 0.1
    @test 0.87 ≤ round(recall; digits=2)
    @test 0.56 ≤ round(precision; digits=2)
    @test 0.69 ≤ round(F_score; digits=2)

    (; labeled_fraction, recall, precision, F_score) = run_and_validate_segmentation(
        first(filter(c -> (c.case_number == 14 && c.satellite == "aqua"), dataset)),
        LopezAcosta2019Tiling.Segment();
    )
    @test 0.33 ≈ labeled_fraction atol = 0.1
    @test 0.85 ≤ round(recall; digits=2)
    @test 0.31 ≤ round(precision; digits=2)
    @test 0.46 ≤ round(F_score; digits=2)

    (; labeled_fraction, recall, precision, F_score) = run_and_validate_segmentation(
        first(filter(c -> (c.case_number == 61 && c.satellite == "aqua"), dataset)),
        LopezAcosta2019Tiling.Segment();
    )
    @test 0.27 ≈ labeled_fraction atol = 0.1
    @test 0.71 ≤ round(recall; digits=2)
    @test 0.67 ≤ round(precision; digits=2)
    @test 0.70 ≤ round(F_score; digits=2)

    (; labeled_fraction, recall, precision, F_score) = run_and_validate_segmentation(
        first(filter(c -> (c.case_number == 63 && c.satellite == "aqua"), dataset)),
        LopezAcosta2019Tiling.Segment();
    )
    @test 0.58 ≈ labeled_fraction atol = 0.1
    @test 0.85 ≤ round(recall; digits=2)
    @test 0.62 ≤ round(precision; digits=2)
    @test 0.73 ≤ round(F_score; digits=2)
end
# Segmentation fault with case 102 with Julia 1.12, passed with Julia 1.11
@testitem "LopezAcosta2019Tiling.Segment - sample of cases" setup = [Segmentation] tags = [
    :e2e
] begin
    using StatsBase: mean
    dataset = Watkins2026Dataset(; ref="v0.1")

    results = run_and_validate_segmentation(
        filter(case -> (case.case_number % 17 == 0), dataset),
        LopezAcosta2019Tiling.Segment();
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
    @test round(mean_F_score; digits=1) ≥ 0.1

    # return current performance
    @show mean_recall
    @show mean_precision
    @show mean_F_score
end

@testitem "LopezAcosta2019Tiling.Segment - image types" tags = [:e2e] setup = [Segmentation] begin
    using Images: RGB, RGBA, n0f8, n6f10, n4f12, n2f14, n0f16, float32, float64
    dataset = Watkins2026Dataset(; ref="v0.1")

    case = first(filter(c -> (c.case_number == 6 && c.satellite == "terra"), dataset))
    algorithm = LopezAcosta2019Tiling.Segment()
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
end
