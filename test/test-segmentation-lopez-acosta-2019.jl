include("./segmentation_utils.jl")
using Images: segment_labels, segment_mean, labels_map

@ntestset "$(@__FILE__)" begin

    # Symbols returned in `_, intermediate_results = LopezAcosta2019()(...; return_intermediate_results=true)`
    # which can be written to PNGs using `save()`
    intermediate_result_image_names = [
        :truecolor,
        :falsecolor,
        :landmask_dilated,
        :landmask_non_dilated,
        :cloudmask,
        :sharpened_truecolor_image,
        :sharpened_gray_truecolor_image,
        :normalized_image,
        :segA,
        :watersheds_segB_product,
        :segF,
        :segment_mean_truecolor,
        :segment_mean_falsecolor,
    ]

    @ntestset "Lopez-Acosta 2019" begin
        @ntestset "Validated data" begin
            data_loader = Watkins2025GitHub(;
                ref="a451cd5e62a10309a9640fbbe6b32a236fcebc70"
            )
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
                result_images_to_save=intermediate_result_image_names,
            )
            @test all(filter(!broken_cases, results).success)
            @test any(filter(broken_cases, results).success) broken = true
        end
        truecolor = load(
            "./test_inputs/pipeline/input_pipeline/20220914.aqua.truecolor.250m.tiff"
        )
        falsecolor = load(
            "./test_inputs/pipeline/input_pipeline/20220914.aqua.reflectance.250m.tiff"
        )
        landmask = load("./test_inputs/pipeline/input_pipeline/landmask.tiff")
        @ntestset "Image types" begin
            region = (200:400, 500:700)
            supported_types = [n0f8, n6f10, n4f12, n2f14, n0f16, float32, float64]
            for target_type in supported_types
                @info "Image type: $target_type"
                intermediate_results_callback = save_results_callback(
                    "./test_outputs/segmentation-LopezAcosta2019-$(target_type)-$(Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS"))";
                    names=intermediate_result_image_names,
                )
                segments = LopezAcosta2019()(
                    target_type.(truecolor[region...]),
                    target_type.(falsecolor[region...]),
                    target_type.(landmask[region...]);
                    intermediate_results_callback,
                )
                @show segments
                @test length(segment_labels(segments)) == 10
            end
        end
        @ntestset "Medium size" begin
            intermediate_results_callback = save_results_callback(
                "./test_outputs/segmentation-LopezAcosta2019-medium-size-$( Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS"))";
                names=intermediate_result_image_names,
            )
            segments = LopezAcosta2019()(
                truecolor, falsecolor, landmask; intermediate_results_callback
            )
            @show segments
            @test length(segment_labels(segments)) == 44
        end
    end
end
