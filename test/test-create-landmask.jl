@testset "Create Landmask" begin
    println("------------------------------------------------")
    println("------------ Create Landmask Test --------------")

    # define constants, maybe move to test config file
    matlab_landmask_file = "$(test_data_dir)/matlab_landmask.png"
    strel_file = "$(test_data_dir)/se.csv"
    struct_elem = readdlm(strel_file, ',', Bool) # read in original matlab structuring element -  a disk-shaped kernel with radius of 50 px
    matlab_landmask = float64.(load(matlab_landmask_file)[lm_test_region...])
    lm_image = float64.(load(landmask_file)[lm_test_region...])
    test_image = float64.(load(truecolor_test_image_file)[lm_test_region...])
    @time landmask = IceFloeTracker.create_landmask(lm_image, struct_elem)

    @time masked_image = IceFloeTracker.apply_landmask(test_image, landmask)

    matlab_landmask_filename =
        "$(test_output_dir)/matlab_landmask_test_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist matlab_landmask matlab_landmask_filename

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
    @test test_similarity(landmask, convert(BitMatrix, matlab_landmask), 0.005)
    #@test (@test_approx_eq_sigma_eps masked_image masked_matlab_image) #TODO #Matlab output is not dilated
end
