@testitem "Normalize Image" begin
    using Images:
        channelview,
        colorview,
        RGB,
        @test_approx_eq_sigma_eps,
        test_approx_eq_sigma_eps,
        strel_diamond,
        float64,
        load
    using Dates: Dates

    include("config.jl")

    struct_elem2 = strel_diamond((5, 5)) #original matlab structuring element -  a disk-shaped kernel with radius of 2 px
    matlab_normalized_img_file = "$(test_data_dir)/matlab_normalized.png"
    matlab_sharpened_file = "$(test_data_dir)/matlab_sharpened.png"
    matlab_diffused_file = "$(test_data_dir)/matlab_diffused.png"
    matlab_equalized_file = "$(test_data_dir)/matlab_equalized.png"
    # flip ocean mask to land mask
    landmask_bitmatrix = convert(
        BitMatrix, float64.(load(current_landmask_file)[test_region...])
    )
    landmask_no_dilate = convert(
        BitMatrix, float64.(load(landmask_no_dilate_file)[test_region...])
    )
    input_image = float64.(load(truecolor_test_image_file)[test_region...])
    matlab_norm_image = float64.(load(matlab_normalized_img_file)[test_region...])
    matlab_sharpened = float64.(load(matlab_sharpened_file))
    matlab_diffused = float64.(load(matlab_diffused_file)[test_region...])
    matlab_equalized = float64.(load(matlab_equalized_file))

    @info "Process Image - Diffusion"
    input_landmasked = apply_landmask(input_image, landmask_no_dilate)

    pmd = PeronaMalikDiffusion(0.1, 0.1, 5, "exponential")
    @time image_diffused = nonlinear_diffusion(input_landmasked, pmd)

    @test (@test_approx_eq_sigma_eps image_diffused matlab_diffused [0, 0] 0.0054) ===
        nothing

    # dmw: not sure what these are for
    @test (@test_approx_eq_sigma_eps input_landmasked image_diffused [0, 0] 0.004) ===
        nothing
    @test (@test_approx_eq_sigma_eps input_landmasked matlab_diffused [0, 0] 0.007) ===
        nothing

    diffused_image_filename =
        "$(test_output_dir)/diffused_test_image_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    @persist image_diffused diffused_image_filename

    @info "Process Image - Equalization"

    ## Equalization
    masked_view = (channelview(matlab_diffused))
    eq = [
        LopezAcosta2019._adjust_histogram(masked_view[i, :, :], 255, 10, 10, 0.86) for
        i in 1:3
    ]
    image_equalized = colorview(RGB, eq...)
    @test (@test_approx_eq_sigma_eps image_equalized matlab_equalized [0, 0] 0.051) ===
        nothing

    equalized_image_filename =
        "$(test_output_dir)/equalized_test_image_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    @persist image_equalized equalized_image_filename

    @info "Process Image - Sharpening"

    ## Sharpening
    @time sharpenedimg = LopezAcosta2019.imsharpen(input_image, landmask_no_dilate)
    @time image_sharpened_gray = LopezAcosta2019.imsharpen_gray(
        sharpenedimg, landmask_bitmatrix
    )
    @test (@test_approx_eq_sigma_eps image_sharpened_gray matlab_sharpened [0, 0] 0.046) ===
        nothing

    sharpened_image_filename =
        "$(test_output_dir)/sharpened_test_image_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    @persist image_sharpened_gray sharpened_image_filename

    @info "Process Image - Normalization"

    ## Normalization
    @time normalized_image = LopezAcosta2019.normalize_image(
        sharpenedimg, image_sharpened_gray, landmask_bitmatrix, struct_elem2
    )

    #test for percent difference in normalized images
    eps = test_approx_eq_sigma_eps(normalized_image, matlab_norm_image, ones(2), 0.1, true)
    @info "Epsilon: " * string(eps)

    @test test_approx_eq_sigma_eps(
        normalized_image, matlab_norm_image, ones(2), 0.1, true
    ) <= 0.05
    nothing

    normalized_image_filename =
        "$(test_output_dir)/normalized_test_image_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    @persist normalized_image normalized_image_filename
end
