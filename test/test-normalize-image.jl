@testset "Normalize Image" begin
    println("-------------------------------------------------")
    println("---------- Create Normalization Test ------------")
    strel_file = "$(test_data_dir)/se.csv"
    struct_elem = readdlm(strel_file, ',', Bool)
    strel_h, strel_w = ceil.(Int, size(struct_elem))
    strel_file2 = "$(test_data_dir)/se2.csv"
    struct_elem2 = readdlm(strel_file2, ',', Bool)
    test_region = (1:2707, 1:4458)
    matlab_normalized_img_file = "$(test_data_dir)/matlab_normalized.tiff"
    landmask = load("$(test_data_dir)/current_landmask.png")
    landmask_bitmatrix = convert(BitMatrix, landmask)
    input_image = load(test_image_file)[test_region...]
    matlab_norm_img = load(matlab_normalized_img_file)[test_region...]
    println("-------------- Process Image ----------------")
    @time normalized_image = IceFloeTracker.normalize_image(
        input_image, landmask_bitmatrix, struct_elem2; kappa=90, clip=0.95
    )

    # test for percent difference in normalized images
    @test (@test_approx_eq_sigma_eps normalized_image[
        strel_h:(end - strel_h), strel_w:(end - strel_w)
    ] matlab_norm_img[strel_h:(end - strel_h), strel_w:(end - strel_w)] [0, 0] 0.058) ==
        nothing
end
