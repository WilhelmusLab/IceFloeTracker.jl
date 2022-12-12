
@testset "Segmentation-A" begin
    println("------------------------------------------------")
    println("------------ Create Segmentation-A Test --------------")

    ice_water_discriminated_image = float64.(load(ice_water_discrim_test_file))
    cloudmask = convert(BitMatrix, load(cloudmask_test_file))
    ice_labels = DelimitedFiles.readdlm("$(test_data_dir)/ice_labels_julia.csv", ',')
    ice_labels = Int64.(vec(ice_labels))
    matlab_segmented_A = float64.(load("$(test_data_dir)/matlab_segmented_A.png"))
    matlab_segmented_A_bitmatrix = convert(BitMatrix, matlab_segmented_A)
    matlab_segmented_ice = convert(
        BitMatrix, float64.(load("$(test_data_dir)/matlab_segmented_ice.png"))
    )

    println("---------- Segment Image - Direct Method ------------")
    @time segmented_ice_filled, segmented_A = IceFloeTracker.segmentation_A(
        ice_water_discriminated_image, cloudmask, ice_labels
    )

    segmented_a_filename =
        "$(test_output_dir)/segmented_a_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist segmented_A segmented_a_filename

    @time segmented_ice = IceFloeTracker.kmeans_segmentation(
        ice_water_discriminated_image, ice_labels
    )

    segmented_ice_filename =
        "$(test_output_dir)/segmented_ice_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist segmented_ice segmented_ice_filename

    @test typeof(segmented_A) == typeof(matlab_segmented_A_bitmatrix)
    @test test_similarity(matlab_segmented_A_bitmatrix, segmented_A, 0.077)

    @test typeof(segmented_ice) == typeof(matlab_segmented_ice)
    @test test_similarity(matlab_segmented_ice, segmented_ice, 0.077)
end
