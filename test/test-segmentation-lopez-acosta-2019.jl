using Images: segment_labels, segment_mean, labels_map

function segment_comparison(;
    name::Union{AbstractString,Nothing}=nothing,
    validated::Union{SegmentedImage,Nothing}=nothing,
    measured::Union{SegmentedImage,Nothing}=nothing,
)
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

    return (;
        name,
        success=true,
        measured_area,
        validated_area,
        fractional_intersection,
        error=nothing,
    )
end

@ntestset "$(@__FILE__)" begin
    @ntestset "Lopez-Acosta 2019" begin
        @ntestset "Validated data" begin
            data_loader = Watkins2025GitHub(;
                ref="a451cd5e62a10309a9640fbbe6b32a236fcebc70"
            )
            @ntestset "visible floes, no clouds, no artifacts" begin
                happy_path_dataset = data_loader(;
                    case_filter=c -> (
                        c.visible_floes == "yes" &&
                        c.cloud_category_manual == "none" &&
                        c.artifacts == "no"
                    ),
                )
                results = []
                for validation_data in happy_path_dataset.data
                    name = validation_data.name
                    validated = validation_data.validated_labeled_floes
                    try
                        measured = LopezAcosta2019()(
                            RGB.(validation_data.modis_truecolor),
                            RGB.(validation_data.modis_falsecolor),
                            validation_data.modis_landmask,
                        )
                        @show measured
                        datestamp = Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS")
                        save(
                            "./test_outputs/segmentation-LopezAcosta2019-$(name)-$(datestamp)-mean-labels.png",
                            map(i -> segment_mean(measured, i), labels_map(measured)),
                        )
                        save(
                            "./test_outputs/segmentation-LopezAcosta2019-$(name)-$(datestamp)-validated-mean-labels.png",
                            map(i -> segment_mean(validated, i), labels_map(validated)),
                        )
                        save(
                            "./test_outputs/segmentation-LopezAcosta2019-$(name)-$(datestamp)-truecolor.png",
                            validation_data.modis_truecolor,
                        )
                        push!(results, segment_comparison(; name, validated, measured))
                    catch e
                        @warn "$(name) failed, $(e)"
                        push!(
                            results,
                            segment_comparison(;
                                name,
                                validated=validation_data.validated_labeled_floes,
                                measured=nothing,
                            ),
                        )
                    end
                end
                results_df = DataFrame(results)
                @info sort(results_df)
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
