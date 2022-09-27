@testset "Segmentation-B" begin
    println("------------------------------------------------")
    println("------------ Create Segmentation-B Test --------------")

    sharpened_image = float64.(load(sharpened_test_image_file))
    segmented_a_ice_mask = convert(BitMatrix, load(segmented_a_ice_mask_file))
    cloudmask = convert(BitMatrix, load(cloudmask_test_file))
    struct_elem2 = readdlm(strel_file_2, ',', Bool)
    matlab_segmented_B = convert(BitMatrix, load("$(test_data_dir)/matlab_segmented_b.png"))

    @time segmented_B = IceFloeTracker.segmentation_B(
        sharpened_image, cloudmask, segmented_a_ice_mask, struct_elem2
    )

    segmented_b_filename =
        "$(test_output_dir)/segmented_b_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist segmented_B segmented_b_filename

    @test typeof(segmented_B) == typeof(matlab_segmented_B)
    @test test_similarity(matlab_segmented_B, segmented_B, 0.078)
end
