@testset "Discriminate Ice-Water" begin
    println("------------------------------------------------")
    println("------------ Create Discrimination Test --------------")

    matlab_Z2 = "$(test_data_dir)/matlab_Z3.png"
    strel_file2 = "$(test_data_dir)/se2.csv"
    struct_elem2 = readdlm(strel_file2, ',', Bool)

    @time landmask = IceFloeTracker.add_padding(
        load(current_landmask_file)[test_region...], Pad(:reflect, (50, 50))
    )
    @time landmask_bitmatrix = convert(BitMatrix, landmask)
    input_image = IceFloeTracker.add_padding(
        load(truecolor_test_image_file)[test_region...], Pad(:reflect, (50, 50))
    )

    @time normalized_image = IceFloeTracker.normalize_image(
        input_image, landmask_bitmatrix, struct_elem2; kappa=90, clip=0.95
    )
    normalized_image = IceFloeTracker.remove_padding(
        normalized_image, Pad((50, 50), (50, 50))
    )

    Z3 = IceFloeTracker.discriminate_ice_water(
        reflectance_test_image_file, normalized_image
    )

    @test (@test_approx_eq_sigma_eps Z3 matlab_Z3 [0, 0] 0.005) == nothing
end
