@testitem "LopezAcosta2019.Segment – simple case" tags = [:e2e, :smoke] begin
    import DataFrames: DataFrame, nrow
    import Images: RGB
    dataset = Watkins2026Dataset(; ref="v0.1")

    case = first(filter(c -> (c.case_number == 6 && c.satellite == "terra"), dataset))
    segments = LopezAcosta2019.Segment()(
        RGB.(modis_truecolor(case)),
        RGB.(modis_falsecolor(case)),
        RGB.(modis_landmask(case)),
    )
    @show segments
    expected_segment_count = validated_floe_properties(case) |> DataFrame |> nrow
    @test length(segments.segment_labels) ≈ expected_segment_count rtol = 0.7
end

@testitem "LopezAcosta2019.Segment – sample of cases" setup = [Segmentation] tags = [:e2e] begin
    dataset = Watkins2026Dataset(; ref="v0.1")
    passing = c -> c.case_number % 11 == 0
    # Case 4 has only very small floes, while case 39 is missing data from Aqua.
    formerly_broken = c -> (c.case_number == 4 || (c.case_number == 39 && c.satellite == "aqua"))
    broken = c -> false  # `broken_cases` once fixed, for regression testing
    results = run_and_validate_segmentation(
        filter(c -> (passing(c) || formerly_broken(c) || broken(c)), dataset),
        LopezAcosta2019.Segment();
        output_directory="./test_outputs/",
    )
    @test all(filter(!broken, results).success)
    @test any(filter(broken, results).success) broken = true
end

@testitem "LopezAcosta2019.Segment – detailed tests" setup = [Segmentation] tags = [:e2e] begin
    dataset = Watkins2026Dataset(; ref="v0.1")
    (; labeled_fraction, recall, precision, F_score) = run_and_validate_segmentation(
        first(filter(c -> (c.case_number == 6 && c.satellite == "terra"), dataset)),
        LopezAcosta2019.Segment();
        output_directory="./test_outputs/",
    )
    @test 0.12 ≈ labeled_fraction atol = 0.1
    @test 0.27 ≤ round(recall; digits=2)
    @test 0.75 ≤ round(precision; digits=2)
    @test 0.40 ≤ round(F_score; digits=2)

    (; labeled_fraction, recall, precision, F_score) = run_and_validate_segmentation(
        first(filter(c -> (c.case_number == 14 && c.satellite == "aqua"), dataset)),
        LopezAcosta2019.Segment();
        output_directory="./test_outputs/",
    )
    @test 0.05 ≈ labeled_fraction atol = 0.1
    @test 0.36 ≤ round(recall; digits=2)
    @test 0.5 ≤ round(precision; digits=2)
    @test 0.46 ≤ round(F_score; digits=2)

    (; labeled_fraction, recall, precision, F_score) = run_and_validate_segmentation(
        first(filter(c -> (c.case_number == 61 && c.satellite == "aqua"), dataset)),
        LopezAcosta2019.Segment();
        output_directory="./test_outputs/",
    )
    @test 0.13 ≈ labeled_fraction atol = 0.1
    @test 0.37 ≤ round(recall; digits=2)
    @test 0.74 ≤ round(precision; digits=2)
    @test 0.50 ≤ round(F_score; digits=2)

    (; labeled_fraction, recall, precision, F_score) = run_and_validate_segmentation(
        first(filter(c -> (c.case_number == 63 && c.satellite == "aqua"), dataset)),
        LopezAcosta2019.Segment();
        output_directory="./test_outputs/",
    )
    @test labeled_fraction ≈ 0.579 rtol = 0.1 broken = true
    @test 0.90 ≤ round(recall; digits=2) broken = true
    @test 0.62 ≤ round(precision; digits=2)
    @test 0.73 ≤ round(F_score; digits=2) broken = true
end

@testitem "LopezAcosta2019.Segment – image types" setup = [Segmentation] tags = [:e2e] begin
    import Images: RGB, n0f8, n6f10, n4f12, n2f14, n0f16, float32, float64
    dataset = Watkins2026Dataset(; ref="v0.1")
    case::Case = first(filter(c -> (c.case_number == 6 && c.satellite == "terra"), dataset))
    algorithm = LopezAcosta2019.Segment()
    baseline = run_and_validate_segmentation(
        case, algorithm; output_directory="./test_outputs/"
    )

    @test results_invariant_for(RGB; baseline, algorithm, case)
    @test results_invariant_for(RGBA; baseline, algorithm, case) broken = true
    @test results_invariant_for(n0f8; baseline, algorithm, case) broken = true
    @test results_invariant_for(n6f10; baseline, algorithm, case) broken = true
    @test results_invariant_for(n4f12; baseline, algorithm, case) broken = true
    @test results_invariant_for(n2f14; baseline, algorithm, case) broken = true
    @test results_invariant_for(n0f16; baseline, algorithm, case) broken = true
    @test results_invariant_for(float32; baseline, algorithm, case) broken = true
    @test results_invariant_for(float64; baseline, algorithm, case) broken = true
    @test results_invariant_for(RGB, n0f8; baseline, algorithm, case)
    @test results_invariant_for(RGB, n6f10; baseline, algorithm, case) broken = true
    @test results_invariant_for(RGB, n4f12; baseline, algorithm, case) broken = true
    @test results_invariant_for(RGB, n2f14; baseline, algorithm, case) broken = true
    @test results_invariant_for(RGB, n0f16; baseline, algorithm, case)
    @test results_invariant_for(RGB, float32; baseline, algorithm, case)
    @test results_invariant_for(RGB, float64; baseline, algorithm, case)
end
