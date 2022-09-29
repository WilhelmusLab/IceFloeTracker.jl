@testset "Segmentation-B" begin
    println("------------------------------------------------")
    println("------------ Create Segmentation-B Test --------------")

    sharpened_image = float64.(load(sharpened_test_image_file))
    segmented_a_ice_mask = convert(BitMatrix, load(segmented_a_ice_mask_file))
    cloudmask = convert(BitMatrix, load(cloudmask_test_file))
    struct_elem2 = readdlm(strel_file_2, ',', Bool)
    matlab_segmented_B_filled = convert(
        BitMatrix, load("$(test_data_dir)/matlab_segmented_b_filled.png")
    )
    matlab_segmented_B_ice = convert(
        BitMatrix, load("$(test_data_dir)/matlab_segmented_b_ice.png")
    )

    @time segmented_B_filled, segmented_B_ice = IceFloeTracker.segmentation_B(
        sharpened_image, cloudmask, segmented_a_ice_mask, struct_elem2
    )

    segmented_b_filled_filename =
        "$(test_output_dir)/segmented_b_filled" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist segmented_B_filled segmented_b_filled_filename

    segmented_b_ice_filename =
        "$(test_output_dir)/segmented_b_ice" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist segmented_B_ice segmented_b_ice_filename

    @test typeof(segmented_B_filled) == typeof(matlab_segmented_B_filled)
    @test test_similarity(matlab_segmented_B_filled, segmented_B_filled, 0.05)

    @test typeof(segmented_B_ice) == typeof(matlab_segmented_B_ice)
    @test test_similarity(matlab_segmented_B_ice, segmented_B_ice, 0.078)
end
