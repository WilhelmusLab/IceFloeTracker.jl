
@testitem "Segmentation-A" begin
    import DelimitedFiles: readdlm
    import Images: float64, load
    import IceFloeTracker.Segmentation: kmeans_binarization
    import IceFloeTracker.LopezAcosta2019: IceDetectionLopezAcosta2019

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
    fc_image = float64.(load(falsecolor_test_image_file)[test_region...])


    println("---------- Segment Image - Direct Method ------------")
    fc_image = load("$(test_data_dir)/beaufort-chukchi-seas_falsecolor.2020162.aqua.250m.tiff")[test_region...]
    fc_landmasked = apply_landmask(fc_image, landmask)
    @time segmented_ice_cloudmasked = LopezAcosta2019.segmented_ice_cloudmasking(
        ice_water_discriminated_image, fc_landmasked, cloudmask
    )

    # Set up the ice detection algorithm (can this be imported, instead?)
    band_7_max=Float64(5 / 255)
    band_2_min=Float64(230 / 255)
    band_1_min=Float64(240 / 255)
    band_7_max_relaxed=Float64(10 / 255)
    band_1_min_relaxed=Float64(190 / 255)
    possible_ice_threshold=Float64(75 / 255)

    @time segmented_ice = kmeans_binarization(ice_water_discriminated_image, fc_landmasked;
        ice_labels_algorithm=IceDetectionLopezAcosta2019())


    segmented_ice = kmeans_binarization(
        ice_water_discriminated_image,
        fc_image;
        k=4,
        maxiter=50,
        random_seed=45,
        ice_labels_algorithm=IceDetectionLopezAcosta2019()
        ) 

    # check: are there any regions that are nonzero under the cloudmask, since it was applied in discriminate ice water?
    segmented_ice_cloudmasked = apply_cloudmask(segmented_ice, cloudmask) 
    @time segmented_A = LopezAcosta2019.clean_binary_floes(segmented_ice_cloudmasked)
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
