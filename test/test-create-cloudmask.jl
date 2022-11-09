@testset "Create Cloudmask" begin
    println("-------------------------------------------------")
    println("------------ Create Cloudmask Test --------------")

    # define constants, maybe move to test config file
    matlab_cloudmask_file = "$(test_data_dir)/matlab_cloudmask.tiff"
    println("--------- Create and apply cloudmask --------")
    ref_image = float64.(load(reflectance_test_image_file)[test_region...])

    matlab_cloudmask = float64.(load(matlab_cloudmask_file))
    @time cloudmask, ref_image_b7 = IceFloeTracker.create_cloudmask(ref_image)
    @time masked_image, clouds_channel = IceFloeTracker.apply_cloudmask(
        ref_image, cloudmask
    )
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
    IceFloeTracker.@persist clouds_channel clouds_channel_filename
    ref_image_b7_filename =
        "$(test_output_dir)/ref_image_b7_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist ref_image_b7 ref_image_b7_filename
    # test for percent difference in landmask images
    @test (@test_approx_eq_sigma_eps masked_image matlab_cloudmask [0, 0] 0.005) === nothing
end
