@testset "Discriminate Ice-Water" begin
    println("------------------------------------------------")
    println("------------ Create Ice-Water Discrimination Test --------------")

    reflectance_image = load(reflectance_test_image_file)[test_region...]
    reflectance_image_band7 = load(reflectance_b7_test_file)
    landmask = load(current_landmask_file)
    landmask_bitmatrix = convert(BitMatrix, landmask)
    normalized_image = load(normalized_test_file)
    clouds_channel = load(clouds_channel_test_file)
    matlab_ice_water_discrim = load("$(test_data_dir)/matlab_ice_water_discrim.png")

    ice_water_discrim = IceFloeTracker.discriminate_ice_water(
        reflectance_image,
        reflectance_image_band7,
        normalized_image,
        landmask_bitmatrix,
        clouds_channel,
    )
    ice_water_discrim_filename =
        "$(test_output_dir)/ice_water_discrim_test_image_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist ice_water_discrim ice_water_discrim_filename

    @test (@test_approx_eq_sigma_eps ice_water_discrim matlab_ice_water_discrim [0, 0] 0.061) ==
        nothing
end
