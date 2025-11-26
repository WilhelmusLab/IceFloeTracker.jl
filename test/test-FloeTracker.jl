@testitem "FloeTracker" begin
    import Dates: DateTime
    import DataFrames: DataFrame
    dataset = Watkins2026Dataset(; ref="v0.1")

    @testset "Basic functionality" begin
        filter!(c -> c.case_number == 6, dataset)
        sort!([:pass_time], dataset)
        segmenter = LopezAcosta2019Tiling.Segment()
        segmentation_results =
            segmenter.(
                modis_truecolor.(dataset),
                modis_falsecolor.(dataset),
                modis_landmask.(dataset),
            )
        tracker = FloeTracker(;
            filter_function=FilterFunction(),
            matching_function=MinimumWeightMatchingFunction(),
        )
        tracking_results = tracker(segmentation_results, info(dataset).pass_time)
        @test isa(tracking_results, DataFrame)
        @test "trajectory_uuid" in names(tracking_results)
    end
end
