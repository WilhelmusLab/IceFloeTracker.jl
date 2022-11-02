
@testset "Segmentation-A" begin
    println("------------------------------------------------")
    println("------------ Create Segmentation-A Test --------------")

    reflectance_image = float64.(load(reflectance_test_image_file)[test_region...])
    landmask = convert(BitMatrix, load(current_landmask_file))
    cloudmask = convert(BitMatrix, load(cloudmask_test_file))
    ice_water_discriminated_image = float64.(load(ice_water_discrim_test_file))
    matlab_segmented_A = float64.(load("$(test_data_dir)/matlab_segmented_A.png"))
    matlab_segmented_A_bitmatrix = convert(BitMatrix, matlab_segmented_A)

    println("---------- Segment Image - Direct Method ------------")
    @time segmented_ice, segmented_ice_filled, segmented_A = IceFloeTracker.segmentation_A(
        reflectance_image, ice_water_discriminated_image, landmask, cloudmask
    )

    segmented_a_filename =
        "$(test_output_dir)/segmented_a_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist segmented_A segmented_a_filename

    segmented_ice_filename =
        "$(test_output_dir)/segmented_ice_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist segmented_ice segmented_ice_filename

    @test typeof(segmented_A) == typeof(matlab_segmented_A_bitmatrix)
    @test test_similarity(matlab_segmented_A_bitmatrix, segmented_A, 0.0845)
    # TODO(tjd) @test_reference "matlab_segmented_A" "segmented_A" by=psnr_equality(20)
end
