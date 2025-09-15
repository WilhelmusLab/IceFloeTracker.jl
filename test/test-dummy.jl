@testitem "dummy functions" begin
    # TODO: remove once dummy functions are deleted
    @test IceFloeTracker.Filtering.dummy_filtering_function()
    @test IceFloeTracker.Morphology.dummy_morphology_function()
    @test IceFloeTracker.Preprocessing.dummy_preprocessing_function()
    @test IceFloeTracker.Segmentation.dummy_segmentation_function()
    @test IceFloeTracker.Tracking.dummy_tracking_function()
    @test IceFloeTracker.Utils.dummy_utils_function()
end