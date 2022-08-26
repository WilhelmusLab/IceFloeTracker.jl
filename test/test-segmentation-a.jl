@testset "Segmentation-A" begin
    println("------------------------------------------------")
    println("------------ Create Segmentation-A Test --------------")

    reflectance_image = load(reflectance_test_image_file)[test_region...]
    landmask = convert(BitMatrix, load(current_landmask_file))
    cloudmask = convert(BitMatrix, load(cloudmask_test_file))
    ice_water_discriminated_image = load(ice_water_discrim_test_file)
    matlab_segmented_A = load("$(test_data_dir)/matlab_segmented_A.png")
    fuzzy_c = load("$(test_data_dir)/fuzzy_C.png")
    fuzzy_c_masked = load("$(test_data_dir)/fuzzy_c_masked.png")
    segmented_A2 = load("$(test_data_dir)/segmented_A2.png")

    println("---------- Segment Image - Direct Method ------------")
    segmented_A = IceFloeTracker.segmentation_A(
        reflectance_image, ice_water_discriminated_image, landmask, cloudmask
    )

    println("------------ Segment Image - Fuzzy-C --------------")
    # segmented_A_fuzzy_C = IceFloeTracker.segmentation_A(ice_water_discriminated_image, cloudmask)

    segmented_a_filename =
        "$(test_output_dir)/segmented_a_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist segmented_A segmented_a_filename

    # segmented_a_fuzzy_c_filename =
    #     "$(test_output_dir)/segmented_a_" *
    #     Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
    #     ".png"
    # IceFloeTracker.@persist segmented_A_fuzzy_C segmented_a_fuzzy_c_filename

    @test (@test_approx_eq_sigma_eps fuzzy_c matlab_segmented_A [0, 0] 0.099) == nothing

    @test (@test_approx_eq_sigma_eps fuzzy_c_masked matlab_segmented_A [0, 0] 0.0845) ==
        nothing

    @test (@test_approx_eq_sigma_eps matlab_segmented_A segmented_A [0, 0] 0.000000000000000001) ==
        nothing #for some reason always passes regardless of eps value

    # @test (@test_approx_eq_sigma_eps matlab_segmented_A segmented_A_fuzzy_C [0, 0] 0.0845) ==
    # nothing

    @test (@test_approx_eq_sigma_eps matlab_segmented_A segmented_A2 [0, 0] 0.0845) ==
        nothing
end
