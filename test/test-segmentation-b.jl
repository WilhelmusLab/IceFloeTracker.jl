@testset "Segmentation-B" begin
    println("------------------------------------------------")
    println("------------ Create Segmentation-B Test --------------")

    sharpened_image = float64.(load(sharpened_test_image_file))
    segmented_a_ice_mask = convert(BitMatrix, load(segmented_a_ice_mask_file))
    cloudmask = convert(BitMatrix, load(cloudmask_test_file))
    matlab_segmented_B_filled = convert(
        BitMatrix, load("$(test_data_dir)/matlab_segmented_b_filled.png")
    )
    matlab_segmented_B_ice = convert(
        BitMatrix, load("$(test_data_dir)/matlab_segmented_b_ice.png")
    )

    matlab_not_ice_mask = float64.(load("$(test_data_dir)/matlab_not_ice_mask.png")) .> 0.5

    @time not_ice_mask, segmented_B_filled, segmented_B_ice = IceFloeTracker.segmentation_B(
        sharpened_image, cloudmask, segmented_a_ice_mask, strel_diamond((3, 3))
    )

    IceFloeTracker.@persist not_ice_mask "./test_outputs/not_ice_mask.png" true

    IceFloeTracker.@persist segmented_B_filled "./test_outputs/segmented_b_filled.png" true

    IceFloeTracker.@persist segmented_B_ice "./test_outputs/segmented_b_ice.png" true

    @test typeof(segmented_B_filled) == typeof(matlab_segmented_B_filled)
    @test test_similarity(matlab_segmented_B_filled, segmented_B_filled, 0.04)

    @test typeof(segmented_B_ice) == typeof(matlab_segmented_B_ice)
    @test test_similarity(matlab_segmented_B_ice, segmented_B_ice, 0.005)

    @test typeof(not_ice_mask) == typeof(matlab_not_ice_mask)
    @test test_similarity(not_ice_mask, (matlab_not_ice_mask), 0.033)
end
