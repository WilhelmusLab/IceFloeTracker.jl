@testset "Segmentation-B" begin
    println("------------------------------------------------")
    println("------------ Create Segmentation-B Test --------------")

    sharpened_image = float64.(load(sharpened_test_image_file))
    segmented_a_ice_mask = convert(BitMatrix, load(segmented_a_ice_mask_file))
    cloudmask = convert(BitMatrix, load(cloudmask_test_file))
    matlab_watershed_intersect = convert(
        BitMatrix, load("$(test_data_dir)/matlab_watershed_intersect.png")
    )
    matlab_not_ice_mask = float64.(load("$(test_data_dir)/matlab_I.png"))
    matlab_not_ice_bit = matlab_not_ice_mask .> 0.499

    @time segB = IceFloeTracker.segmentation_B(
        sharpened_image, cloudmask, segmented_a_ice_mask
    )

    IceFloeTracker.@persist segB.not_ice "./test_outputs/segB_not_ice_mask.png" true
    IceFloeTracker.@persist segB.watershed_intersect "./test_outputs/segB_watershed.png" true
    IceFloeTracker.@persist matlab_not_ice_mask "./test_outputs/matlab_not_ice_mask.png" true
    IceFloeTracker.@persist matlab_watershed_intersect "./test_outputs/matlab_watershed.png" true

    @test typeof(segB.not_ice) == typeof(matlab_not_ice_mask)
    @test (@test_approx_eq_sigma_eps segB.not_ice matlab_not_ice_mask [0, 0] 0.001) ==
        nothing

    @test test_similarity(matlab_watershed_intersect, segB.watershed_intersect, 0.011)
end
