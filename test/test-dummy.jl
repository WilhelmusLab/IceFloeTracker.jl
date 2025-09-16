@testitem "dummy functions" begin
    # TODO: remove once dummy functions are deleted
    @test !isempty(IceFloeTracker.Filtering.dummy_filtering_function())
    @test !isempty(IceFloeTracker.Morphology.dummy_morphology_function())
    @test !isempty(IceFloeTracker.Preprocessing.dummy_preprocessing_function())
    @test !isempty(IceFloeTracker.Segmentation.dummy_segmentation_function())
    @test !isempty(IceFloeTracker.Tracking.dummy_tracking_function())
    @test !isempty(IceFloeTracker.Utils.dummy_utils_function())
end