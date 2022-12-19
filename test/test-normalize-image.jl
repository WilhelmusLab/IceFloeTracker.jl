@testset "Normalize Image" begin
    println("-------------------------------------------------")
    println("---------- Create Normalization Test ------------")

    matlab_normalized_img_file = "$(test_data_dir)/matlab_normalized.png"
    matlab_sharpened_file = "$(test_data_dir)/matlab_sharpened.png"
    matlab_diffused_file = "$(test_data_dir)/matlab_diffused.png"
    # matlab_gammared_file = "$(test_data_dir)/matlab_gammared.png"
    # matlab_gammagreen_file = "$(test_data_dir)/matlab_gammagreen.png"
    matlab_equalized_file = "$(test_data_dir)/matlab_equalized.png"
    landmask_bitmatrix = convert(BitMatrix, float64.(load(current_landmask_file)))
    landmask_no_dilate = convert(BitMatrix, float64.(load(landmask_no_dilate_file)))
    input_image = float64.(load(truecolor_test_image_file)[test_region...])
    input_image = IceFloeTracker.apply_landmask(input_image, landmask_no_dilate)
    matlab_norm_image = float64.(load(matlab_normalized_img_file)[test_region...])
    matlab_sharpened = float64.(load(matlab_sharpened_file))
    matlab_diffused = float64.(load(matlab_diffused_file)[test_region...])
    # matlab_gammared = float64.(load(matlab_gammared_file)[ice_floe_test_region...])
    # matlab_gammagreen = float64.(load(matlab_gammagreen_file)[ice_floe_test_region...])
    matlab_equalized = float64.(load(matlab_equalized_file))

    println("-------------- Process Image - Diffusion ----------------")

    ## Diffusion
    @time image_diffused = IceFloeTracker.diffusion(input_image, 0.1, 75, 3)

    @test (@test_approx_eq_sigma_eps image_diffused matlab_diffused [0, 0] 0.0054) ==
        nothing

    @test (@test_approx_eq_sigma_eps input_image image_diffused [0, 0] 0.004) == nothing
    @test (@test_approx_eq_sigma_eps input_image matlab_diffused [0, 0] 0.007) == nothing

    diffused_image_filename =
        "$(test_output_dir)/diffused_test_image_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist image_diffused diffused_image_filename

    println("-------------- Process Image - Equalization ----------------")

    ## Equalization
    masked_view = (channelview(image_diffused))
    eq = [
        IceFloeTracker._adjust_histogram(masked_view[i, :, :], 255, 22, 22, 0.87) for
        i in 1:3
    ]
    image_equalized = colorview(RGB, eq...)
    @test (@test_approx_eq_sigma_eps image_equalized matlab_equalized [0, 0] 0.056) ==
        nothing

    equalized_image_filename =
        "$(test_output_dir)/equalized_test_image_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist image_equalized equalized_image_filename

    println("-------------- Process Image - Sharpening ----------------")

    ## Sharpening
    @time sharpenedimg = IceFloeTracker.imsharpen(input_image, landmask_no_dilate)
    @time image_sharpened_gray = IceFloeTracker.imsharpen_gray(
        sharpenedimg, landmask_bitmatrix
    )
    @test (@test_approx_eq_sigma_eps image_sharpened_gray matlab_sharpened [0, 0] 0.052) ==
        nothing

    sharpened_image_filename =
        "$(test_output_dir)/sharpened_test_image_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist image_sharpened_gray sharpened_image_filename

    println("-------------- Process Image - Normalization ----------------")

    ## Normalization
    @time normalized_image = IceFloeTracker.normalize_image(
        sharpenedimg, image_sharpened_gray, landmask_bitmatrix
    )

    #test for percent difference in normalized images
    @test (@test_approx_eq_sigma_eps normalized_image matlab_norm_image [0, 0] 0.056) ==
        nothing

    normalized_image_filename =
        "$(test_output_dir)/normalized_test_image_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist normalized_image normalized_image_filename
end
