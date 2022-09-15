@testset "Normalize Image" begin
    println("-------------------------------------------------")
    println("---------- Create Normalization Test ------------")
    strel_file2 = "$(test_data_dir)/se2.csv"
    struct_elem2 = readdlm(strel_file2, ',', Bool)
    matlab_normalized_img_file = "$(test_data_dir)/matlab_normalized.tiff"
    landmask = float32.(load(current_landmask_file))
    landmask_bitmatrix = convert(BitMatrix, landmask)
    input_image = float32.(load(truecolor_test_image_file)[test_region...])
    matlab_norm_image = float32.(load(matlab_normalized_img_file)[test_region...])
    println("-------------- Process Image ----------------")
    @time normalized_image = IceFloeTracker.normalize_image(
        input_image, landmask_bitmatrix, struct_elem2; kappa=90, clip=0.95
    )
    normalized_image_filename =
        "$(test_output_dir)/normalized_test_image_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist normalized_image normalized_image_filename

    # test for percent difference in normalized images
    @test (@test_approx_eq_sigma_eps normalized_image matlab_norm_image [0, 0] 0.058) ==
        nothing
end
