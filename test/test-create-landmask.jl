@testset "Create Landmask" begin
    println("------------------------------------------------")
    println("------------ Create Landmask Test --------------")

    # define constants, maybe move to test config file
    matlab_landmask_file = "$(test_data_dir)/matlab_landmask.png"
    strel_file = "$(test_data_dir)/se.csv"
    lm_test_region = (3000:3750, 1000:1550)
    num_pixels_closing = 50
    struct_elem = readdlm(strel_file, ',', Bool)
    strel_h, strel_w = ceil.(Int, size(struct_elem) ./ 2)
    lm_image = IceFloeTracker.add_padding(
        load(landmask_file)[lm_test_region...], Pad(:replicate, (50, 50))
    )
    matlab_landmask = load(matlab_landmask_file)[lm_test_region...]
    test_image = IceFloeTracker.add_padding(
        load(truecolor_test_image_file)[lm_test_region...], Pad(:replicate, (50, 50))
    )

    @time landmask = IceFloeTracker.create_landmask(
        lm_image, struct_elem; num_pixels_closing=num_pixels_closing
    )
    println("--------- Run second time for JIT warmup --------")
    @time landmask = IceFloeTracker.create_landmask(
        lm_image, struct_elem; num_pixels_closing=num_pixels_closing
    )

    @time masked_image = IceFloeTracker.apply_landmask(test_image, landmask)

    landmask_filename =
        "$(test_output_dir)/landmask_test_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    landmask = IceFloeTracker.@persist IceFloeTracker.remove_padding(
        landmask, Pad((50, 50), (50, 50))
    ) landmask_filename

    masked_image_filename =
        "$(test_output_dir)/landmasked_truecolor_test_image_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    masked_image = IceFloeTracker.@persist IceFloeTracker.remove_padding(
        masked_image, Pad((50, 50), (50, 50))
    ) masked_image_filename

    # test for percent difference in landmask images, ignore edges because we are not padding in Julia before applying strel_file
    @test (@test_approx_eq_sigma_eps landmask matlab_landmask [0, 0] 0.005) == nothing

    # TO DO: add test of applied landmask
end
