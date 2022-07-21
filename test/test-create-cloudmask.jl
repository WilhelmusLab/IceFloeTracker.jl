@testset "Create Cloudmask" begin
    println("-------------------------------------------------")
    println("------------ Create Cloudmask Test --------------")

    # define constants, maybe move to test config file
    matlab_cloudmask_file = "$(test_data_dir)/matlab_cloudmask.tiff"
    println("--------- Create and apply cloudmask --------")
    ref_image = IceFloeTracker.add_padding(
        load(reflectance_test_image_file)[test_region...], Pad(:replicate, (50, 50))
    )
    matlab_cloudmask = load(matlab_cloudmask_file)
    @time cloudmask = IceFloeTracker.create_cloudmask(ref_image)
    @time masked_image, clouds_channel = IceFloeTracker.apply_cloudmask(
        ref_image, cloudmask
    )
    masked_image = IceFloeTracker.remove_padding(masked_image, Pad((50, 50), (50, 50)))

    # test for percent difference in landmask images
    @test (@test_approx_eq_sigma_eps masked_image matlab_cloudmask [0, 0] 0.005) == nothing
end
