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
    @time landmask_no_dilate = landmask.non_dilated

    @time masked_image = IceFloeTracker.apply_landmask(test_image, landmask.dilated)
    @time masked_image_no_dilate = IceFloeTracker.apply_landmask(
        test_image, .!landmask_no_dilate
    )

    # test for percent difference in landmask images
    @test test_similarity(.!landmask.dilated, convert(BitMatrix, matlab_landmask), 0.005)
    @test test_similarity(
        landmask.non_dilated, convert(BitMatrix, matlab_landmask_no_dilate), 0.005
    )

    # test for in-place allocation reduction
    @time normal_lm = IceFloeTracker.apply_landmask(test_image, landmask.dilated)
    @time IceFloeTracker.apply_landmask!(test_image, landmask.dilated)

    x = @allocated IceFloeTracker.apply_landmask(test_image, landmask.dilated)
    @info("normal allocated: $x")
    y = @allocated IceFloeTracker.apply_landmask!(test_image, landmask.dilated)
    @info("in-place allocated: $y")
    @test x > y

    # test that the test image has been updated in-place and equals the new image with landmask applied
    @test(test_image == normal_lm)

    # persist imgs
    matlab_landmask_filename =
        "$(test_output_dir)/matlab_landmask_test_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist matlab_landmask matlab_landmask_filename

    landmask_filename =
        "$(test_output_dir)/landmask_test_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist landmask.dilated landmask_filename

    landmask_no_dilate_filename =
        "$(test_output_dir)/landmask_test_no_dilate_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist landmask.non_dilated landmask_no_dilate_filename

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
end
