@testset "Normalize Image" begin
    println("-------------------------------------------------")
    println("---------- Create Normalization Test ------------")
    struct_elem2 = strel_diamond((5,5)) #original matlab structuring element -  a disk-shaped kernel with radius of 2 px
    matlab_normalized_img_file = "$(test_data_dir)/matlab_normalized.tiff"
    matlab_sharpened_file = "$(test_data_dir)/matlab_sharpened.png"
    landmask = float64.(load(current_landmask_file))
    landmask_bitmatrix = convert(BitMatrix, landmask)
    input_image = float64.(load(truecolor_test_image_file)[test_region...])
    matlab_norm_image = float64.(load(matlab_normalized_img_file)[test_region...])
    matlab_sharpened = load(matlab_sharpened_file)[ice_floe_test_region...]

    println("-------------- Process Image ----------------")
    @time sharpenedimg = IceFloeTracker.imsharpen(input_image)
    @time image_sharpened_gray = IceFloeTracker.imsharpen_gray(sharpenedimg, landmask_bitmatrix)
    @time normalized_image = IceFloeTracker.normalize_image(sharpenedimg, image_sharpened_gray, landmask_bitmatrix, struct_elem2)

    # test method with default se
    @test IceFloeTracker.normalize_image(sharpenedimg, image_sharpened_gray, landmask_bitmatrix) == normalized_image
    
    # test for percent difference in normalized images
    @test (@test_approx_eq_sigma_eps normalized_image matlab_norm_image [0, 0] 0.058) ==
        nothing

    @test (@test_approx_eq_sigma_eps image_sharpened_gray[ice_floe_test_region...] matlab_sharpened [
        0, 0
    ] 0.065) == nothing

    normalized_image_filename =
        "$(test_output_dir)/normalized_test_image_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist normalized_image normalized_image_filename

    sharpened_image_filename =
        "$(test_output_dir)/sharpened_test_image_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist image_sharpened_gray sharpened_image_filename
end
