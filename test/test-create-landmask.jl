@testset "Create Landmask" begin
    println("------------------------------------------------")
    println("------------ Create Landmask Test --------------")

    # define constants, maybe move to test config file
    matlab_landmask_file = "$(test_data_dir)/matlab_landmask.png"
    matlab_landmask_no_dilate_file = "$(test_data_dir)/matlab_landmask_no_dilate.png"
    strel_file = "$(test_data_dir)/se.csv"
    struct_elem = readdlm(strel_file, ',', Bool) # read in original matlab structuring element -  a disk-shaped kernel with radius of 50 px
    matlab_landmask = float64.(load(matlab_landmask_file)[lm_test_region...])
    matlab_landmask_no_dilate =
        float64.(load(matlab_landmask_no_dilate_file)[lm_test_region...])
    lm_image = float64.(load(landmask_file)[lm_test_region...])
    test_image = load(truecolor_test_image_file)[lm_test_region...]
    @time landmask = IceFloeTracker.create_landmask(lm_image, struct_elem)

    # Test method with default se
    @test landmask == IceFloeTracker.create_landmask(lm_image)

    # Generate testing files
    @time landmask_no_dilate = IceFloeTracker.create_landmask_bool(lm_image)

    @time masked_image = IceFloeTracker.apply_landmask(test_image, landmask)
    @time masked_image_no_dilate = IceFloeTracker.apply_landmask_no_dilate(
        test_image, .!landmask_no_dilate
    )

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

    landmask_no_dilate_filename =
        "$(test_output_dir)/landmask_test_no_dilate_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist landmask_no_dilate landmask_no_dilate_filename

    masked_image_filename =
        "$(test_output_dir)/landmasked_truecolor_test_image_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist masked_image masked_image_filename

    masked_image_no_dilate_filename =
        "$(test_output_dir)/landmasked_truecolor_test_no_dilate_image_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist masked_image_no_dilate masked_image_no_dilate_filename

    # test for percent difference in landmask images
    @test test_similarity(.!landmask, convert(BitMatrix, matlab_landmask), 0.005)
    @test test_similarity(
        landmask_no_dilate, convert(BitMatrix, matlab_landmask_no_dilate), 0.005
    )
end
