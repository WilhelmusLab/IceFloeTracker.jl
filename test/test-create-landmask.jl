@testset "Create Landmask" begin
    println("------------------------------------------------")
    println("------------ Create Landmask Test --------------")

    # define constants, maybe move to test config file
    matlab_landmask_file = "$(test_data_dir)/matlab_landmask.png"
    strel_file = "$(test_data_dir)/se.csv"
    struct_elem = readdlm(strel_file, ',', Bool)
    matlab_landmask = float32.(load(matlab_landmask_file)[lm_test_region...])
    lm_image = float32.(load(landmask_file)[lm_test_region...])
    test_image = float32.(load(truecolor_test_image_file)[lm_test_region...])
    @time landmask = IceFloeTracker.create_landmask(lm_image, struct_elem)

    @time masked_image = IceFloeTracker.apply_landmask(test_image, landmask)

    landmask_filename =
        "$(test_output_dir)/landmask_test_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist landmask landmask_filename

    masked_image_filename =
        "$(test_output_dir)/landmasked_truecolor_test_image_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist masked_image masked_image_filename

    # test for percent difference in landmask images
    @test (@test_approx_eq_sigma_eps landmask matlab_landmask [0, 0] 0.005) == nothing

    #@test (@test_approx_eq_sigma_eps masked_image masked_matlab_image) #TODO #Matlab output is not dilated
end
