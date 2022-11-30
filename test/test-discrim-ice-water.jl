@testset "Discriminate Ice-Water" begin
    println("------------------------------------------------")
    println("------------ Create Ice-Water Discrimination Test --------------")
    input_image = float64.(load(truecolor_test_image_file)[test_region...])
    reflectance_image = float64.(load(reflectance_test_image_file)[test_region...])
    landmask =  convert(BitMatrix, load(current_landmask_file))
    cloudmask = IceFloeTracker.create_cloudmask(reflectance_image)
    matlab_ice_water_discrim =
        float64.(load("$(test_data_dir)/matlab_ice_water_discrim.png"))

    image_sharpened = IceFloeTracker.imsharpen(input_image)
    ice_water_discrim = 
    IceFloeTracker.discriminate_ice_water(
        reflectance_image, 
        image_sharpened, landmask, cloudmask)
    @test (@test_approx_eq_sigma_eps ice_water_discrim matlab_ice_water_discrim [0, 0] 0.061) ==
    nothing

    # persist generated image
    ice_water_discrim_filename =
        "$(test_output_dir)/ice_water_discrim_test_image_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist ice_water_discrim ice_water_discrim_filename
end
