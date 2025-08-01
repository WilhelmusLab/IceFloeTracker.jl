using Images: segment_labels, segment_mean, labels_map

function segmentation_comparison(;
    validated::Union{SegmentedImage,Nothing}=nothing,
    measured::Union{SegmentedImage,Nothing}=nothing,
)::NamedTuple
    if !isnothing(validated)
        validated_binary = Gray.(labels_map(validated) .> 0)
        validated_area = sum(channelview(validated_binary))
    else
        validated_area = missing
    end

    if !isnothing(measured)
        measured_binary = Gray.(labels_map(measured) .> 0)
        measured_area = sum(channelview(measured_binary))
    else
        measured_area = missing
    end

    if !isnothing(validated) && !isnothing(measured)
        intersection = Gray.(channelview(measured_binary) .&& channelview(validated_binary))
        fractional_intersection =
            fractional_intersection = sum(channelview(intersection)) / validated_area
    else
        fractional_intersection = missing
    end

    return (; measured_area, validated_area, fractional_intersection)
end

function run_segmentation_over_multiple_cases(
    cases,
    algorithm::IceFloeSegmentationAlgorithm;
    output_path::AbstractString="./test_outputs",
    save_outputs::Bool=true,
)
    results = []
    for case::ValidationDataCase in cases
        validated = case.validated_labeled_floes
        name = case.name
        let measured, success, error
            @info "starting $(name)"
            try
                measured = algorithm(
                    RGB.(case.modis_truecolor),
                    RGB.(case.modis_falsecolor),
                    case.modis_landmask,
                )
                @info "$(name) succeeded"
                success = true
                error = nothing
            catch error
                @warn "$(name) failed: $(error)"
                success = false
                measured = nothing
            end

            # Store the aggregate results
            push!(
                results,
                merge(
                    (; name, success, error), segmentation_comparison(; validated, measured)
                ),
            )

            if save_outputs
                datestamp = Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS")
                !isnothing(measured) && save(
                    joinpath(
                        output_path,
                        "segmentation-$(typeof(algorithm))-$(name)-$(datestamp)-mean-labels.png",
                    ),
                    map(i -> segment_mean(measured, i), labels_map(measured)),
                )
                !isnothing(validated) && save(
                    joinpath(
                        output_path,
                        "segmentation-$(typeof(algorithm))-$(name)-$(datestamp)-validated-mean-labels.png",
                    ),
                    map(i -> segment_mean(validated, i), labels_map(validated)),
                )
                !isnothing(case.modis_truecolor) && save(
                    joinpath(
                        output_path,
                        "segmentation-$(typeof(algorithm))-$(name)-$(datestamp)-truecolor.png",
                    ),
                    case.modis_truecolor,
                )
            end
        end
    end
    results_df = DataFrame(results)
    return results_df
end

@ntestset "$(@__FILE__)" begin
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
                    dataset.data, LopezAcosta2019()
                )
                @info results

                # Run tests on aggregate results
                # First be sure we have the right number of results
                # if this fails, then everything below needs re-evaluating
                @test nrow(results) == nrow(dataset.metadata)

                # Check the number of successes
                rows_successful = subset(results, :success => ByRow(==(true)))
                @test nrow(rows_successful) >= 6
                nrow(rows_successful) > 6 && @warn "new passing cases"
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
                    dataset.data, LopezAcosta2019()
                )
                @info results

                # Run tests on aggregate results
                # First be sure we have the right number of results
                @test nrow(results) == nrow(dataset.metadata)

                # Check the number of successes
                rows_successful = subset(results, :success => ByRow(==(true)))
                @test nrow(rows_successful) >= 6
                nrow(rows_successful) > 6 && @warn "new passing cases"
            end
        end

        truecolor = load(
            "./test_inputs/pipeline/input_pipeline/20220914.aqua.truecolor.250m.tiff"
        )
        falsecolor = load(
            "./test_inputs/pipeline/input_pipeline/20220914.aqua.reflectance.250m.tiff"
        )
        landmask = load("./test_inputs/pipeline/input_pipeline/landmask.tiff")
        @ntestset "Smoke test" begin
            region = (200:400, 500:700)
            supported_types = [n0f8, n6f10, n4f12, n2f14, n0f16, float32, float64]
            for target_type in supported_types
                @info "Image type: $target_type"
                segments = LopezAcosta2019()(
                    target_type.(truecolor[region...]),
                    target_type.(falsecolor[region...]),
                    target_type.(landmask[region...]),
                )
                @show segments
                save(
                    "./test_outputs/segmentation-Lopez-Acosta-2019-mean-labels_$(target_type)_$(Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS")).png",
                    map(i -> segment_mean(segments, i), labels_map(segments)),
                )
                @test length(segment_labels(segments)) == 10
            end
        end

        @ntestset "Full size" begin
            segments = LopezAcosta2019()(truecolor, falsecolor, landmask)
            @show segments
            save(
                "./test_outputs/segmentation-Lopez-Acosta-2019-mean-labels_$(Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS")).png",
                map(i -> segment_mean(segments, i), labels_map(segments)),
            )
            @test length(segment_labels(segments)) == 44
        end
    end
end
