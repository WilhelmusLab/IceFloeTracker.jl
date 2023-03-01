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

    @time segB = IceFloeTracker.segmentation_B(
        sharpened_image, cloudmask, segmented_a_ice_mask
    )

    IceFloeTracker.@persist segB.not_ice "./test_outputs/segB_not_ice_mask.png" true

    IceFloeTracker.@persist segB.filled "./test_outputs/segB_filled.png" true

    IceFloeTracker.@persist segB.ice "./test_outputs/segB_ice.png" true

    @test typeof(segB.filled) == typeof(matlab_segmented_B_filled)
    @test test_similarity(matlab_segmented_B_filled, segB.filled, 0.04)

    @test typeof(segB.ice) == typeof(matlab_segmented_B_ice)
    @test test_similarity(matlab_segmented_B_ice, segB.ice, 0.005)

    @test typeof(segB.not_ice) == typeof(matlab_not_ice_mask)
    @test test_similarity(segB.not_ice, matlab_not_ice_mask, 0.033)
end
