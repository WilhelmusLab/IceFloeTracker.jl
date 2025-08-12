include("./segmentation_utils.jl")
using Images: segment_labels, segment_mean, labels_map

@ntestset "$(@__FILE__)" begin
    @ntestset "Lopez-Acosta 2019" begin
        data_loader = Watkins2025GitHub(; ref="a451cd5e62a10309a9640fbbe6b32a236fcebc70")
        @ntestset "Sample of cases" begin
            broken_cases =
                c -> (c.case_number == 4 || (c.case_number == 39 && c.satellite == "aqua"))
            formerly_broken_cases = c -> false  # `broken_cases` once fixed, for regression testing
            results = run_segmentation_over_multiple_cases(
                data_loader,
                c -> (
                    c.case_number % 17 == 0 || formerly_broken_cases(c) || broken_cases(c)
                ),
                LopezAcosta2019();
                output_directory="./test_outputs/",
            )
            @test all(filter(!broken_cases, results).success)
            @test any(filter(broken_cases, results).success) broken = true
        end
        @ntestset "Detailed tests" begin
            check_case(
                LopezAcosta2019(),
                data_loader,
                c -> (c.case_number == 6 && c.satellite == "terra");
                labeled_fraction_goal=0.119,
                recall_goal=0.315,
                precision_goal=0.770,
                F_score_goal=0.447,
            )
            check_case(
                LopezAcosta2019(),
                data_loader,
                c -> (c.case_number == 14 && c.satellite == "aqua");
                labeled_fraction_goal=0.052,
                recall_goal=0.360,
                precision_goal=0.857,
                F_score_goal=0.507,
            )
            check_case(
                LopezAcosta2019(),
                data_loader,
                c -> (c.case_number == 61 && c.satellite == "aqua");
                labeled_fraction_goal=0.132,
                recall_goal=0.379,
                precision_goal=0.754,
                F_score_goal=0.504,
            )
            (; labeled_fraction, recall, precision, F_score) = run_segmentation_over_one_case(
                data_loader,
                c -> (c.case_number == 63 && c.satellite == "aqua"),
                LopezAcosta2019(),
            )
            @test labeled_fraction ≈ 0.579 rtol = 0.1 broken = true
            @test 0.901 ≤ recall broken = true
            @test 0.620 ≤ precision broken = true
            @test 0.734 ≤ F_score broken = true
        end
        @ntestset "Image types" begin
            case = first(
                data_loader(;
                    case_filter=c -> (c.case_number == 6 && c.satellite == "terra")
                ),
            )
            supported_types = [n0f8, n6f10, n4f12, n2f14, n0f16, float32, float64]
            for target_type in supported_types
                @info "Image type: $target_type"
                intermediate_results_callback = save_results_callback(
                    "./test_outputs/segmentation-LopezAcosta2019-$(target_type)-$(Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS"))";
                )
                segments = LopezAcosta2019()(
                    target_type.(RGB.(case.modis_truecolor)),
                    target_type.(RGB.(case.modis_falsecolor)),
                    target_type.(RGB.(case.modis_landmask));
                    intermediate_results_callback,
                )
                @show segments
                @test length(segment_labels(segments)) ≈ 65 atol = 2
            end
        end
    end
end
