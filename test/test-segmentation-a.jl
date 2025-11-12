
@testitem "Segmentation-A" begin
    import DelimitedFiles: readdlm
    import Images: float64, load

    include("config.jl")
    include("test_error_rate.jl")

    ice_water_discriminated_image =
        float64.(load("$(test_data_dir)/matlab_ice_water_discrim.png"))
    cloudmask = .!convert(BitMatrix, load(cloudmask_test_file))
    landmask = convert(BitMatrix, load(current_landmask_file)[test_region...]) # Test landmask file has ocean == 1
    ice_labels = readdlm("$(test_data_dir)/ice_labels_julia.csv", ',')
    ice_labels = Int64.(vec(ice_labels))
    matlab_segmented_A = float64.(load("$(test_data_dir)/matlab_segmented_A.png"))
    matlab_segmented_A_bitmatrix = convert(BitMatrix, matlab_segmented_A)
    matlab_segmented_ice = convert(
        BitMatrix, float64.(load("$(test_data_dir)/matlab_segmented_ice.png"))
    )
    matlab_segmented_ice_cloudmasked =
        load("$(test_data_dir)/matlab_segmented_ice_cloudmasked.png") .> 0.5

    println("---------- Segment Image - Direct Method ------------")
    @time segmented_ice_cloudmasked = LopezAcosta2019.segmented_ice_cloudmasking(
        ice_water_discriminated_image, cloudmask, ice_labels
    )

    @time segmented_A = LopezAcosta2019.segmentation_A(segmented_ice_cloudmasked)

    @time segmented_ice = kmeans_segmentation(ice_water_discriminated_image, ice_labels)

    @test typeof(segmented_A) == typeof(matlab_segmented_A_bitmatrix)
    @test test_similarity(matlab_segmented_A_bitmatrix, segmented_A, 0.039)

    @test typeof(segmented_ice_cloudmasked) == typeof(matlab_segmented_ice_cloudmasked)
    @test test_similarity(
        convert(BitMatrix, apply_landmask(segmented_ice_cloudmasked, landmask)),
        matlab_segmented_ice_cloudmasked,
        0.051,
    )

    @test typeof(segmented_ice) == typeof(matlab_segmented_ice)
    @test test_similarity(matlab_segmented_ice, segmented_ice, 0.001)
end
