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

    # Persist output images
    cloudmask_filename =
        "$(test_output_dir)/cloudmask_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist cloudmask cloudmask_filename
    masked_image_filename =
        "$(test_output_dir)/cloudmasked_reflectance_test_image_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist masked_image masked_image_filename
    clouds_channel_filename =
        "$(test_output_dir)/clouds_channel_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist clds_channel clouds_channel_filename
end
