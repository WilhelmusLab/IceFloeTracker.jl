@testset "Create Cloudmask" begin
    println("-------------------------------------------------")
    println("------------ Create Cloudmask Test --------------")

    # define constants, maybe move to test config file
    ref_image_file = "$(test_data_dir)/cloudmask_test_image.tiff"
    matlab_cloudmask_file = "$(test_data_dir)/matlab_cloudmask.tiff"
    println("--------- Create and apply cloudmask --------")
    ref_image = load(ref_image_file)
    matlab_cloudmask = load(matlab_cloudmask_file)
    @time cloudmask = IceFloeTracker.create_cloudmask(ref_image)
    @time masked_image = IceFloeTracker.apply_cloudmask(ref_image, cloudmask)

    # test for percent difference in landmask images
    @test (@test_approx_eq_sigma_eps masked_image matlab_cloudmask [0, 0] 0.005) == nothing
end
