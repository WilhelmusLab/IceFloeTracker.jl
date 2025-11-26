@testitem "FloeTracker" begin
    import Dates: DateTime
    dataset = Watkins2026Dataset(; ref="v0.1")

    @testset "Basic functionality" begin
        cases = filter(c -> c.case_number == 6, dataset)
        segmenter = LopezAcosta2019Tiling.Segment()
        segmentation_results =
            segmenter.(
                modis_truecolor.(cases), modis_falsecolor.(cases), modis_landmask.(cases)
            )
        @info segmentation_results

        tracker = FloeTracker(;
            filter_function=FilterFunction(),
            matching_function=MinimumWeightMatchingFunction(),
        )
        tracking_results = tracker(segmentation_results, info(cases).pass_time)
        @info tracking_results
    end
end
