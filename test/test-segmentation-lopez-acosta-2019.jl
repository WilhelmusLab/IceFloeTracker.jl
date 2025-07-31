using Images: segment_labels, segment_mean, labels_map

@ntestset "$(@__FILE__)" begin
    @ntestset "Lopez-Acosta 2019" begin
        @ntestset "Validated data" begin
            data_loader = Watkins2025GitHub(;
                cache_dir="./__temp__/Watkins2025GitHub",
                ref="a451cd5e62a10309a9640fbbe6b32a236fcebc70",
            )
            dataset = data_loader()
            results = []
            for validation_data in dataset.data
                name = validation_data.name
                try
                    segments = LopezAcosta2019()(
                        RGB.(validation_data.modis_truecolor),
                        RGB.(validation_data.modis_falsecolor),
                        validation_data.modis_landmask,
                    )
                    @show segments
                    save(
                        "./test_outputs/segmentation-LopezAcosta2019-$(name)_$(Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS"))-mean-labels.png",
                        map(i -> segment_mean(segments, i), labels_map(segments)),
                    )
                    if !isnothing(validation_data.validated_labeled_floes)
                        area = labels_map(segments) .> 0
                        validated_area =
                            labels_map(validation_data.validated_labeled_floes) .> 0
                    else
                        area = nothing
                        validated_area = nothing
                    end

                    push!(
                        results, (; name, success=true, area, validated_area, error=nothing)
                    )
                catch e
                    @warn "$(name) failed, $(e)"
                    push!(results, (; name, success=false, error=e))
                end
                results_df = DataFrame(results)
                @info sort(results_df)
            end
        end

        # truecolor = load(
        #     "./test_inputs/pipeline/input_pipeline/20220914.aqua.truecolor.250m.tiff"
        # )
        # falsecolor = load(
        #     "./test_inputs/pipeline/input_pipeline/20220914.aqua.reflectance.250m.tiff"
        # )
        # landmask = load("./test_inputs/pipeline/input_pipeline/landmask.tiff")
        # @ntestset "Smoke test" begin
        #     region = (200:400, 500:700)
        #     supported_types = [n0f8, n6f10, n4f12, n2f14, n0f16, float32, float64]
        #     for target_type in supported_types
        #         @info "Image type: $target_type"
        #         segments = LopezAcosta2019()(
        #             target_type.(truecolor[region...]),
        #             target_type.(falsecolor[region...]),
        #             target_type.(landmask[region...]),
        #         )
        #         @show segments
        #         save(
        #             "./test_outputs/segmentation-Lopez-Acosta-2019-mean-labels_$(target_type)_$(Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS")).png",
        #             map(i -> segment_mean(segments, i), labels_map(segments)),
        #         )
        #         @test length(segment_labels(segments)) == 10
        #     end
        # end

        # @ntestset "Full size" begin
        #     segments = LopezAcosta2019()(truecolor, falsecolor, landmask)
        #     @show segments
        #     save(
        #         "./test_outputs/segmentation-Lopez-Acosta-2019-mean-labels_$(Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS")).png",
        #         map(i -> segment_mean(segments, i), labels_map(segments)),
        #     )
        #     @test length(segment_labels(segments)) == 44
        # end
    end
end
