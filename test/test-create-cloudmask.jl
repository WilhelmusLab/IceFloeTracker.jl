@testset "Create Cloudmask" begin
    println("-------------------------------------------------")
    println("------------ Create Cloudmask Test --------------")

    # define constants, maybe move to test config file
    matlab_cloudmask_file = "$(test_data_dir)/matlab_cloudmask.tiff"
    println("--------- Create and apply cloudmask --------")
    ref_image = float64.(load(falsecolor_test_image_file)[test_region...])

    matlab_cloudmask = float64.(load(matlab_cloudmask_file))
    @time cloudmask = IceFloeTracker.create_cloudmask(ref_image)
    @time masked_image = IceFloeTracker.apply_cloudmask(ref_image, cloudmask)

    # test for percent difference in cloudmask images
    @test (@test_approx_eq_sigma_eps masked_image matlab_cloudmask [0, 0] 0.005) === nothing

    # test for create_clouds_channel
    clouds_channel_expected = load(clouds_channel_test_file)
    clds_channel = IceFloeTracker.create_clouds_channel(cloudmask, ref_image)
    @test (@test_approx_eq_sigma_eps (clds_channel) (clouds_channel_expected) [0, 0] 0.005) ===
        nothing

    @info "Test image that loads as RGBA"
    pth_RGBA_tiff = "$(test_data_dir)/466-sea_of_okhostk-100km-20040421.terra.truecolor.250m.tiff"
    ref_image = load(pth_RGBA_tiff)
    @test typeof(ref_image) <: Matrix{RGBA{N0f8}}
    cloudmask = IceFloeTracker.create_cloudmask(ref_image)
    @test sum(cloudmask) === 0 # all pixels are clouds
end
