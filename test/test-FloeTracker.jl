@testitem "FloeTracker" begin
    import Dates: DateTime
    data_loader = Watkins2026Dataset(; ref="v0.1")

    @testset "Basic functionality" begin
        cases = data_loader(c -> c.case_number == 6)
        segmenter = LopezAcosta2019Tiling.Segment()
        segmentation_results = [
            segmenter(c.modis_truecolor, c.modis_falsecolor, c.modis_landmask) for
            c in cases
        ]
        @info segmentation_results

        tracker = FloeTracker(;
            filter_function=FilterFunction(),
            matching_function=MinimumWeightMatchingFunction(),
        )
        tracking_results = tracker(segmentation_results, info(cases).pass_time)
        @info tracking_results
    end
end
