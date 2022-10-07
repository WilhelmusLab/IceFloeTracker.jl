@testset "Segmentation-D,E" begin
    println("------------------------------------------------")
    println("------------ Create Segmentation-D,E Test --------------")

    segmented_c = convert(BitMatrix, load(segmented_c_test_file)) ## C3 in matlab
    not_ice_mask = convert(BitMatrix, load(not_ice_mask_test_file))

    watershed_B, watershed_C = IceFloeTracker.segmentation_D_E(not_ice_mask, segmented_c)

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

    #@test typeof(segmented_C) == typeof(matlab_segmented_C)
    #@test test_similarity(matlab_segmented_C, segmented_C, 0.078)
end
