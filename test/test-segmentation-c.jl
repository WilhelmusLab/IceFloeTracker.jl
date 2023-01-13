@testset "Segmentation-C" begin
    println("------------------------------------------------")
    println("------------ Create Segmentation-C Test --------------")

    matlab_segmented_C = convert(BitMatrix, load("$(test_data_dir)/matlab_segmented_c.png"))
    segmented_B_filled = convert(BitMatrix, load(segmented_b_filled_test_file))
    segmented_B_ice = convert(BitMatrix, load(segmented_b_ice_test_file))

    segmented_C = IceFloeTracker.segmentation_C(segmented_B_filled, segmented_B_ice)

    IceFloeTracker.@persist segmented_C "./test_outputs/segmented_c.png" true

    @test typeof(segmented_C) == typeof(matlab_segmented_C)
    @test test_similarity(matlab_segmented_C, segmented_C, 0.001)
end
