include("./segmentation_utils.jl")
using Images: segment_labels, segment_mean, labels_map

@ntestset "$(@__FILE__)" begin

    # Symbols returned in `_, intermediate_results = LopezAcosta2019()(...; return_intermediate_results=true)`
    # which can be written to PNGs using `save()`
    intermediate_result_image_names = [
        :landmask_dilated,
        :landmask_non_dilated,
        :cloudmask,
        :sharpened_truecolor_image,
        :sharpened_gray_truecolor_image,
        :normalized_image,
        :segA,
        :watersheds_segB_product,
        :segF,
    ]

    @ntestset "Lopez-Acosta 2019" begin
        @ntestset "Validated data" begin
            data_loader = Watkins2025GitHub(;
                ref="a451cd5e62a10309a9640fbbe6b32a236fcebc70"
            )
            @ntestset "visible floes, no clouds, no artifacts" begin
                dataset = data_loader(;
                    case_filter=c -> (
                        c.visible_floes == "yes" &&
                        c.cloud_category_manual == "none" &&
                        c.artifacts == "no"
                    ),
                )
                @info dataset.metadata

                results = run_segmentation_over_multiple_cases(
                    dataset.data, LopezAcosta2019(); output_directory="./test_outputs/"
                )
                @info results

                # Run tests on aggregate results
                # First be sure we have the right number of results
                @test nrow(results) == nrow(dataset.metadata)

                # Now check that all cases run through without crashing
                successes = subset(results, :success => ByRow(==(true)))
                @test nrow(results) == nrow(successes) broken = true

                # If not everything works, at least check that we're not introducing new crashes
                expected_successes = 6
                successes = subset(results, :success => ByRow(==(true)))
                @test nrow(successes) >= expected_successes
                nrow(successes) > expected_successes &&
                    @warn "new passing cases: $(nrow(successes)) (update `expected_successes`)"
            end

            @ntestset "visible floes, thin clouds, no artifacts" begin
                dataset = data_loader(;
                    case_filter=c -> (
                        c.visible_floes == "yes" &&
                        c.cloud_category_manual == "thin" &&
                        c.artifacts == "no" &&
                        c.case_number % 5 == 0
                    ),
                )
                @info dataset.metadata
                results = run_segmentation_over_multiple_cases(
                    dataset.data, LopezAcosta2019(); output_directory="./test_outputs/"
                )
                @info results

                # Run tests on aggregate results
                # First be sure we have the right number of results
                @test nrow(results) == nrow(dataset.metadata)

                # Now check that all cases run through without crashing
                successes = subset(results, :success => ByRow(==(true)))
                @test nrow(results) == nrow(successes) broken = true
                # If not everything works, at least check that we're not introducing new crashes
                expected_successes = 3
                @test nrow(successes) >= expected_successes
                nrow(successes) > expected_successes &&
                    @warn "new passing cases: $(nrow(successes)) (update `expected_successes`)"
            end
            @ntestset "random sample" begin
                dataset = data_loader(; case_filter=c -> (c.case_number % 17 == 0))
                @info dataset.metadata
                results = run_segmentation_over_multiple_cases(
                    dataset.data, LopezAcosta2019(); output_directory="./test_outputs/"
                )
                @info results

                # Run tests on aggregate results
                # First be sure we have the right number of results
                @test nrow(results) == nrow(dataset.metadata)

                # Now check that all cases run through without crashing
                successes = subset(results, :success => ByRow(==(true)))
                @test nrow(results) == nrow(successes) broken = true
                # If not everything works, at least check that we're not introducing new crashes
                expected_successes = 3
                @test nrow(successes) >= expected_successes
                nrow(successes) > expected_successes &&
                    @warn "new passing cases: $(nrow(successes)) (update `expected_successes`)"
            end
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
                segments, intermediate_results = LopezAcosta2019()(
                    target_type.(truecolor[region...]),
                    target_type.(falsecolor[region...]),
                    target_type.(landmask[region...]);
                    return_intermediate_results=true,
                )
                datestamp = Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS")
                save_intermediate_images(
                    "./test_outputs/segmentation-Lopez-Acosta-2019-$(target_type)-$(datestamp)-intermediate-results/",
                    intermediate_results;
                    names=intermediate_result_image_names,
                )
                @show segments
                save(
                    "./test_outputs/segmentation-Lopez-Acosta-2019-$(target_type)-$(datestamp)-mean-labels.png",
                    map(i -> segment_mean(segments, i), labels_map(segments)),
                )
                @test length(segment_labels(segments)) == 10
            end
        end
        @ntestset "Medium size" begin
            segments, intermediate_results = LopezAcosta2019()(
                truecolor, falsecolor, landmask; return_intermediate_results=true
            )
            @show segments
            datestamp = Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS")
            save_intermediate_images(
                "./test_outputs/segmentation-LopezAcosta2019-medium-size-$(datestamp)-intermediate-results/",
                intermediate_results;
                names=intermediate_result_image_names,
            )
            save(
                "./test_outputs/segmentation-LopezAcosta2019-medium-size-$(datestamp)-mean-labels.png",
                map(i -> segment_mean(segments, i), labels_map(segments)),
            )
            @test length(segment_labels(segments)) == 44
        end
    end
end
