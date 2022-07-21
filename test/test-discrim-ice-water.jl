@testset "Discriminate Ice-Water" begin
    println("------------------------------------------------")
    println("------------ Create Discrimination Test --------------")

    matlab_Z = "$(test_data_dir)/matlab_Z.png"
    matlab_Z2 = "$(test_data_dir)/matlab_Z2.png"
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

    Z, Z2 = IceFloeTracker.discriminate_ice_water(normalized_image, ref)

    @test (@test_approx_eq_sigma_eps Z matlab_Z [0, 0] 0.005) == nothing

    @test (@test_approx_eq_sigma_eps Z2 matlab_Z2 [0, 0] 0.005) == nothing
end
