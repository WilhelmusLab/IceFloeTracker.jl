@testitem "Preprocess Image" begin
    using IceFloeTracker:
        PeronaMalikDiffusion,
        apply_landmask!,
        unsharp_mask
    using Images:
        channelview,
        colorview,
        RGB,
        Gray,
        @test_approx_eq_sigma_eps,
        test_approx_eq_sigma_eps,
        strel_diamond,
        float64,
        load,
        adjust_histogram,
        complement,
        mreconstruct,
        AdaptiveEqualization,
        dilate
    using Dates: Dates
    using IceFloeTracker.Filtering: _channelwise_adapthisteq
    using IceFloeTracker: LopezAcosta2019
    

    include("config.jl")

    struct_elem2 = strel_diamond((5, 5)) #original matlab structuring element -  a disk-shaped kernel with radius of 2 px
    matlab_normalized_img_file = "$(test_data_dir)/matlab_normalized.png"
    matlab_sharpened_file = "$(test_data_dir)/matlab_sharpened.png"
    matlab_diffused_file = "$(test_data_dir)/matlab_diffused.png"
    matlab_equalized_file = "$(test_data_dir)/matlab_equalized.png"
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

    falsecolor_image = float64.(load(falsecolor_test_image_file)[test_region...])
    cloudmask = IceFloeTracker.create_cloudmask(falsecolor_image) # reversed cloudmask
    matlab_ice_water_discrim =
        float64.(load("$(test_data_dir)/matlab_ice_water_discrim.png"))

    # diffuse -> equalize -> sharpen -> grayscale
    @info "Process Image - Diffusion"
    input_landmasked = apply_landmask(input_image, landmask_no_dilate)

    pmd = PeronaMalikDiffusion(0.1, 0.1, 5, "exponential")
    @time image_diffused = nonlinear_diffusion(input_landmasked, pmd)

    @test (@test_approx_eq_sigma_eps image_diffused matlab_diffused [0, 0] 0.0054) ===
        nothing

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
    image_equalized = _channelwise_adapthisteq(matlab_diffused)
    @test (@test_approx_eq_sigma_eps image_equalized matlab_equalized [0, 0] 0.051) ===
        nothing

    equalized_image_filename =
        "$(test_output_dir)/equalized_test_image_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    @persist image_equalized equalized_image_filename

    @info "Process Image - Sharpening"

    ## Sharpening
    # The functions here could get wrapped together into a "Preprocessing" functor.
    # The imsharpen section is diffusion -> adapt hist eq -> unsharp mask -> Gray
    apply_landmask!(image_diffused, landmask_no_dilate)

    sharpened_truecolor_image = _channelwise_adapthisteq(image_diffused)

    # unsharp masking
    image_sharpened_gray = unsharp_mask(Gray.(sharpened_truecolor_image), 10, 2)
    apply_landmask!(image_sharpened_gray, landmask_bitmatrix)

    @test (@test_approx_eq_sigma_eps image_sharpened_gray matlab_sharpened [0, 0] 0.046) ===
        nothing

    sharpened_image_filename =
        "$(test_output_dir)/sharpened_test_image_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    @persist image_sharpened_gray sharpened_image_filename

    @info "Grayscale Reconstruction"
    # input may need to be Gray.(sharpened_truecolor_image) in case the landmask dilation matters
    markers = complement.(dilate(image_sharpened_gray, strel_diamond((5, 5))))
    mask = complement.(image_sharpened_gray)
    # reconstruction of the complement: floes are dark, leads are bright
    @time reconst_gray = mreconstruct(dilate, markers, mask, strel_diamond((5, 5)))
    apply_landmask!(reconst_gray, landmask_bitmatrix)

    #test for percent difference in normalized images
    eps = test_approx_eq_sigma_eps(reconst_gray, matlab_norm_image, ones(2), 0.1, true)
    @info "Epsilon: " * string(eps)

    @test test_approx_eq_sigma_eps(reconst_gray, matlab_norm_image, ones(2), 0.1, true) <=
        0.05
    nothing

    @info "Ice-water discrimination"

    ice_water_discrim = LopezAcosta2019.discriminate_ice_water(
        falsecolor_image, reconst_gray, landmask_bitmatrix, cloudmask
    )
    @test (@test_approx_eq_sigma_eps ice_water_discrim matlab_ice_water_discrim [0, 0] 0.065) ===
        nothing

    # Which image should be persisted? Equalized grayscale? Reconstructed? Ice-water discrimination?
    preprocessed_image_filename =
        "$(test_output_dir)/preprocessed_test_image_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    @persist reconst_gray preprocessed_image_filename
end

# TODO: Move discrim ice water here so that we don't waste time re-computing variables.