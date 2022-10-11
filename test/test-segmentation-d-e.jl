@testset "Segmentation-D,E" begin
    println("------------------------------------------------")
    println("------------ Create Segmentation-D,E Test --------------")

    matlab_watershed_B = convert(BitMatrix, load("$(test_data_dir)/matlab_watershed_b.png"))
    matlab_watershed_C = convert(BitMatrix, load("$(test_data_dir)/matlab_watershed_c.png"))
    segmented_c = convert(BitMatrix, load(segmented_c_test_file))
    not_ice_mask = convert(BitMatrix, load(not_ice_mask_test_file))

    @time watershed_B, watershed_C = IceFloeTracker.segmentation_D_E(
        not_ice_mask, segmented_c
    )

    watershed_B_filename =
        "$(test_output_dir)/watershed_b_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist watershed_B watershed_B_filename

    watershed_C_filename =
        "$(test_output_dir)/watershed_c_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist watershed_C watershed_C_filename

    @test typeof(watershed_B) == typeof(matlab_watershed_B)
    @test typeof(watershed_C) == typeof(matlab_watershed_C)
    @test test_similarity(matlab_watershed_B, watershed_B, 0.097)
    @test test_similarity(matlab_watershed_C, watershed_C, 0.13)
end
