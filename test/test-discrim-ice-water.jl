@testitem "Discriminate Ice-Water" begin
    using Dates: format, now
    using Images: @test_approx_eq_sigma_eps, float64, load, colorview, Gray

    include("config.jl")

    input_image = float64.(load(truecolor_test_image_file)[test_region...])
    falsecolor_image = float64.(load(falsecolor_test_image_file)[test_region...])
    # Flip the imported landmasks, since it has ocean=0 (i.e. they are ocean masks).
    landmask = convert(BitMatrix, load(current_landmask_file)[test_region...])
    landmask_no_dilate = convert(
        BitMatrix, float64.(load(landmask_no_dilate_file)[test_region...])
    )
    cloudmask = create_cloudmask(falsecolor_image) # reversed cloudmask
    matlab_ice_water_discrim =
        float64.(load("$(test_data_dir)/matlab_ice_water_discrim.png"))

    image_sharpened = LopezAcosta2019.imsharpen(input_image, landmask_no_dilate)
    image_sharpened_gray::Matrix{Gray{Float64}} = colorview(
        Gray, apply_landmask(image_sharpened, landmask)
    )
    normalized_image = LopezAcosta2019.normalize_image(
        image_sharpened, image_sharpened_gray, landmask
    )
    ice_water_discrim = LopezAcosta2019.discriminate_ice_water(
        falsecolor_image, normalized_image, landmask, cloudmask
    )
    @test (@test_approx_eq_sigma_eps ice_water_discrim matlab_ice_water_discrim [0, 0] 0.065) ===
        nothing

    # persist generated image
    ice_water_discrim_filename =
        "$(test_output_dir)/ice_water_discrim_test_image_" *
        format(now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    @persist ice_water_discrim ice_water_discrim_filename
end
