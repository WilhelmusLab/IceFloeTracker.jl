@testset "Normalize Image" begin
    println("-------------------------------------------------")
    println("---------- Create Normalization Test ------------")
    # strel_file = "$(test_data_dir)/se.csv"
    # struct_elem = readdlm(strel_file, ',', Bool)
    # strel_h, strel_w = ceil.(Int, size(struct_elem))
    strel_file2 = "$(test_data_dir)/se2.csv"
    struct_elem2 = readdlm(strel_file2, ',', Bool)
    matlab_normalized_img_file = "$(test_data_dir)/matlab_normalized.tiff"
    landmask = IceFloeTracker.add_padding(
        load(current_landmask_file)[test_region...], Pad(:reflect, (50, 50))
    )
    landmask_bitmatrix = convert(BitMatrix, landmask)
    input_image = IceFloeTracker.add_padding(
        load(truecolor_test_image_file)[test_region...], Pad(:reflect, (50, 50))
    )
    matlab_norm_img = load(matlab_normalized_img_file)[test_region...]
    println("-------------- Process Image ----------------")
    @time normalized_image = IceFloeTracker.normalize_image(
        input_image, landmask_bitmatrix, struct_elem2; kappa=90, clip=0.95
    )
    normalized_image = IceFloeTracker.remove_padding(
        normalized_image, Pad((50, 50), (50, 50))
    )

    # test for percent difference in normalized images
    @test (@test_approx_eq_sigma_eps normalized_image matlab_norm_img [0, 0] 0.058) ==
        nothing
end
