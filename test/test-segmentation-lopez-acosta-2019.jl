@testitem "LopezAcosta2019 – simple case" tags = [:e2e, :smoke] begin
    data_loader = Watkins2025GitHub(; ref="a451cd5e62a10309a9640fbbe6b32a236fcebc70")

    case = first(data_loader(c -> (c.case_number == 6 && c.satellite == "terra")))
    segments = LopezAcosta2019()(
        RGB.(case.modis_truecolor), RGB.(case.modis_falsecolor), RGB.(case.modis_landmask)
    )
    @show segments
    expected_segment_count = case.validated_floe_properties |> DataFrame |> nrow
    @test length(segments.segment_labels) ≈ expected_segment_count rtol = 0.7
end

@testitem "LopezAcosta2019 – sample of cases" setup = [Segmentation] tags = [:e2e] begin
    data_loader = Watkins2025GitHub(; ref="a451cd5e62a10309a9640fbbe6b32a236fcebc70")
    passing = c -> c.case_number % 17 == 0
    broken = c -> (c.case_number == 4 || (c.case_number == 39 && c.satellite == "aqua"))
    formerly_broken = c -> false  # `broken_cases` once fixed, for regression testing
    results = run_and_validate_segmentation(
        data_loader(c -> (passing(c) || formerly_broken(c) || broken(c))),
        LopezAcosta2019();
        output_directory="./test_outputs/",
    )
    @test all(filter(!broken, results).success)
    @test any(filter(broken, results).success) broken = true
end

@testitem "LopezAcosta2019 – detailed tests" setup = [Segmentation] tags = [:e2e] begin
    data_loader = Watkins2025GitHub(; ref="a451cd5e62a10309a9640fbbe6b32a236fcebc70")
    (; labeled_fraction, recall, precision, F_score) = run_and_validate_segmentation(
        first(data_loader(c -> (c.case_number == 6 && c.satellite == "terra"))),
        LopezAcosta2019();
        output_directory="./test_outputs/",
    )
    @test 0.119 ≈ labeled_fraction atol = 0.1
    @test 0.315 ≤ recall
    @test 0.770 ≤ precision
    @test 0.447 ≤ F_score

    (; labeled_fraction, recall, precision, F_score) = run_and_validate_segmentation(
        first(data_loader(c -> (c.case_number == 14 && c.satellite == "aqua"))),
        LopezAcosta2019();
        output_directory="./test_outputs/",
    )
    @test 0.052 ≈ labeled_fraction atol = 0.1
    @test 0.360 ≤ recall
    @test 0.857 ≤ precision
    @test 0.507 ≤ F_score

    (; labeled_fraction, recall, precision, F_score) = run_and_validate_segmentation(
        first(data_loader(c -> (c.case_number == 61 && c.satellite == "aqua"))),
        LopezAcosta2019();
        output_directory="./test_outputs/",
    )
    @test 0.132 ≈ labeled_fraction atol = 0.1
    @test 0.379 ≤ recall
    @test 0.754 ≤ precision
    @test 0.504 ≤ F_score

    (; labeled_fraction, recall, precision, F_score) = run_and_validate_segmentation(
        first(data_loader(c -> (c.case_number == 63 && c.satellite == "aqua"))),
        LopezAcosta2019();
        output_directory="./test_outputs/",
    )
    @test labeled_fraction ≈ 0.579 rtol = 0.1 broken = true
    @test 0.901 ≤ recall broken = true
    @test 0.620 ≤ precision
    @test 0.734 ≤ F_score broken = true
end

@testitem "LopezAcosta2019 – image types" setup = [Segmentation] tags = [:e2e] begin
    using Images: RGB, n0f8, n6f10, n4f12, n2f14, n0f16, float32, float64
    data_loader = Watkins2025GitHub(; ref="a451cd5e62a10309a9640fbbe6b32a236fcebc70")
    case::ValidationDataCase = first(
        data_loader(c -> (c.case_number == 6 && c.satellite == "terra"))
    )
    algorithm = LopezAcosta2019()
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
