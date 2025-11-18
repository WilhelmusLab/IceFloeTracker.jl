@testitem "FloeTracker" begin
    import Dates: DateTime
    data_loader = Watkins2025GitHub(; ref="a451cd5e62a10309a9640fbbe6b32a236fcebc70")

    @testset "Basic functionality" begin
        cases = data_loader(c -> c.case_number == 6)
        segmenter = LopezAcosta2019Tiling.Segment()
        segmentation_results = [
            segmenter(c.modis_truecolor, c.modis_falsecolor, c.modis_landmask) for
            c in cases
        ]
        @info segmentation_results

        tracker = FloeTracker(;
            filter_function=ChainedFilterFunction(),
            matching_function=MinimumWeightMatchingFunction(),
        )
        tracking_results = tracker(
            segmentation_results,
            DateTime.(cases.metadata.start_date), # TODO: return a DateTime from the data loader
        )
        @info tracking_results
    end
end
