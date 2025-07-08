using Images: segment_labels, segment_mean, labels_map

@ntestset "$(@__FILE__)" begin
    @ntestset "Lopez-Acosta 2019" begin
        loadimg(s::String) = float64.(load(s))

        truecolor = loadimg(
            "./test_inputs/pipeline/input_pipeline/20220914.aqua.truecolor.250m.tiff"
        )

        falsecolor = loadimg(
            "./test_inputs/pipeline/input_pipeline/20220914.aqua.reflectance.250m.tiff"
        )

        landmask = loadimg("./test_inputs/pipeline/input_pipeline/landmask.tiff")

        segments = LopezAcosta2019()(truecolor, falsecolor, landmask)

        @show segments
        IceFloeTracker.@persist(
            map(i -> segment_mean(segments, i), labels_map(segments)),
            "./test_outputs/segmentation-Lopez-Acosta-2019-mean-labels.png",
            true
        )
        @test length(segment_labels(segments)) == 44
    end
end
