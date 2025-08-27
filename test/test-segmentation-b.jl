@testitem "Segmentation-B" setup = [Paths, Similarity] begin
    using IceFloeTracker: @test_approx_eq_sigma_eps

    sharpened_image = float64.(load(sharpened_test_image_file))
    segmented_a_ice_mask = convert(BitMatrix, load(segmented_a_ice_mask_file))
    cloudmask = convert(BitMatrix, load(cloudmask_test_file))
    matlab_ice_intersect = convert(
        BitMatrix, load("$(test_data_dir)/matlab_segmented_c.png")
    )
    matlab_not_ice_mask = float64.(load("$(test_data_dir)/matlab_I.png"))
    #matlab_not_ice_bit = matlab_not_ice_mask .> 0.499

    @time segB = IceFloeTracker.segmentation_B(
        sharpened_image, .!cloudmask, segmented_a_ice_mask
    )

    IceFloeTracker.@persist segB.not_ice "./test_outputs/segB_not_ice_mask.png" true
    IceFloeTracker.@persist segB.ice_intersect "./test_outputs/segB_ice_mask.png" true
    IceFloeTracker.@persist matlab_not_ice_mask "./test_outputs/matlab_not_ice_mask.png" true
    IceFloeTracker.@persist matlab_ice_intersect "./test_outputs/matlab_ice_intersect.png" true

    @test typeof(segB.not_ice) == typeof(matlab_not_ice_mask)
    @test (@test_approx_eq_sigma_eps segB.not_ice matlab_not_ice_mask [0, 0] 0.001) ===
        nothing

    @test typeof(segB.not_ice_bit) == typeof(matlab_not_ice_mask .> 0.499)
    @test test_similarity((matlab_not_ice_mask .> 0.499), segB.not_ice_bit, 0.001)

    @test typeof(segB.ice_intersect) == typeof(matlab_ice_intersect)
    @test test_similarity(matlab_ice_intersect, segB.ice_intersect, 0.005)
end
