@testset "Discriminate Ice-Water" begin
    println("------------------------------------------------")
    println("------------ Create Discrimination Test --------------")

    reflectance_image = load(reflectance_test_image_file)[test_region...]
    reflectance_image_band7 = load(reflectance_b7_test_file)
    landmask = load(current_landmask_file)
    landmask_bitmatrix = convert(BitMatrix, landmask)
    normalized_image = load(normalized_test_file)
    clouds_channel = load(clouds_channel_test_file)
    matlab_Z3 = load("$(test_data_dir)/matlab_Z3.png")

    Z3 = IceFloeTracker.discriminate_ice_water(
        reflectance_image,
        reflectance_image_band7,
        normalized_image,
        landmask_bitmatrix,
        clouds_channel,
    )
    Z3_filename =
        "$(test_output_dir)/Z3_test_image_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist Z3 Z3_filename

    @test (@test_approx_eq_sigma_eps Z3 matlab_Z3 [0, 0] 0.18) == nothing
end
