using IceFloeTracker
using Images
using Test
using DelimitedFiles
using Dates

@testset "IceFloeTracker.jl" begin
    test_data_dir = "./data"
    test_image_file = "$(test_data_dir)/NE_Greenland_truecolor.2020162.aqua.250m.tiff"

    @testset "Create Landmask" begin
        println("------------------------------------------------")
        println("------------ Create Landmask Test --------------")

        # define constants, maybe move to test config file
        landmask_file = "$(test_data_dir)/landmask.tiff"
        matlab_landmask_file = "$(test_data_dir)/matlab_landmask.png"
        strel_file = "$(test_data_dir)/se.csv"
        test_region = (3000:3750, 1000:1550)
        num_pixels_closing = 50
        struct_elem = readdlm(strel_file, ',', Bool)
        strel_h, strel_w = ceil.(Int, size(struct_elem) ./ 2)
        lm_image = load(landmask_file)[test_region...]
        matlab_landmask = load(matlab_landmask_file)[test_region...]
        test_image = load(test_image_file)[test_region...]

        @time landmask = IceFloeTracker.create_landmask(
            lm_image, struct_elem; num_pixels_closing=num_pixels_closing
        )
        println("--------- Run second time for JIT warmup --------")
        @time landmask = IceFloeTracker.create_landmask(
            lm_image, struct_elem; num_pixels_closing=num_pixels_closing
        )

        @time masked_image = IceFloeTracker.apply_landmask(test_image, landmask)

        # test for percent difference in landmask images, ignore edges because we are not padding in Julia before applying strel_file
        @test (@test_approx_eq_sigma_eps landmask[
            strel_h:(end - strel_h), strel_w:(end - strel_w)
        ] matlab_landmask[strel_h:(end - strel_h), strel_w:(end - strel_w)] [0, 0] 0.005) ==
            nothing

        # TO DO: add test of applied landmask
    end

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
        @test (@test_approx_eq_sigma_eps masked_image matlab_cloudmask [0, 0] 0.005) ==
            nothing
    end

    @testset "Normalize Image" begin
        println("-------------------------------------------------")
        println("---------- Create Normalization Test ------------")
        strel_file = "$(test_data_dir)/se.csv"
        struct_elem = readdlm(strel_file, ',', Bool)
        strel_h, strel_w = ceil.(Int, size(struct_elem))
        strel_file2 = "$(test_data_dir)/se2.csv"
        struct_elem2 = readdlm(strel_file2, ',', Bool)
        test_region = (1:2707, 1:4458)
        matlab_normalized_img_file = "$(test_data_dir)/matlab_normalized.tiff"
        landmask = load("$(test_data_dir)/current_landmask.png")
        landmask_bitmatrix = convert(BitMatrix, landmask)
        input_image = load(test_image_file)[test_region...]
        matlab_norm_img = load(matlab_normalized_img_file)[test_region...]
        println("-------------- Process Image ----------------")
        @time normalized_image = IceFloeTracker.normalize_image(
            input_image, landmask_bitmatrix, struct_elem2; kappa=90, clip=0.95
        )

        # test for percent difference in normalized images
        @test (@test_approx_eq_sigma_eps normalized_image[
            strel_h:(end - strel_h), strel_w:(end - strel_w)
        ] matlab_norm_img[strel_h:(end - strel_h), strel_w:(end - strel_w)] [0, 0] 0.058) ==
            nothing
    end

    @testset "persist.jl" begin
        println("-------------------------------------------------")
        println("---------- Persist Image Tests ------------")
        img_path = "/landmask.tiff"
        outimage_path = "outimage1.tiff"
        img = Images.load.(test_data_dir * img_path)

        # Test filename in variable
        IceFloeTracker.@persist img outimage_path
        @test isfile(outimage_path)

        # Test filename as string literal
        IceFloeTracker.@persist img "outimage2.tiff"
        @test isfile("outimage2.tiff")

        # Test no-filename call. Default filename startswith 'persisted_mask-' 
        # First clear all files that start with this prefix, if any
        [rm(f) for f in readdir() if startswith(f, "persisted_mask-")]
        @assert length([f for f in readdir() if startswith(f, "persisted_mask-")]) == 0

        IceFloeTracker.@persist img
        @test length([f for f in readdir() if startswith(f, "persisted_mask-")]) == 1

        # clean up - part 1!
        rm(outimage_path)
        rm("outimage2.tiff")

        # Part 2
        IceFloeTracker.@persist identity(img) outimage_path
        IceFloeTracker.@persist identity(img) "outimage2.tiff"
        IceFloeTracker.@persist identity(img) # no file name given
        @test isfile(outimage_path)
        @test isfile("outimage2.tiff")
        @test length([f for f in readdir() if startswith(f, "persisted_mask-")]) == 2

        # Clean up
        rm(outimage_path)
        rm("outimage2.tiff")
        [rm(f) for f in readdir() if startswith(f, "persisted_mask-")]
    end
end
