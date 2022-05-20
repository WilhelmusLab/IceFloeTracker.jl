using IceFloeTracker
using Images
using Test
using DelimitedFiles

@testset "IceFloeTracker.jl" begin

    @testset "Create Landmask" begin
        println("------------------------------------------------")
        println("------------ Create Landmask Test --------------")

        # define constants, maybe move to test config file
        test_data_dir = "./data"
        landmask_file = """$(test_data_dir)/landmask.tiff"""
        matlab_landmask_file = """$(test_data_dir)/matlab_landmask.png"""
        strel_file = """$(test_data_dir)/se.csv"""
        test_region = (3000:3750, 1000:1550)
        num_pixels_closing = 50

        lm_image = load(landmask_file)[test_region...]
        matlab_landmask = load(matlab_landmask_file)[test_region...]
        struct_elem = readdlm(strel_file, ',', Bool)
        strel_h, strel_w = ceil.(Int, size(struct_elem)./2)
        @time masked_image = IceFloeTracker.create_landmask(lm_image, struct_elem; num_pixels_closing=num_pixels_closing)
        println("----------Run second time for JIT warmup---------")
        @time masked_image = IceFloeTracker.create_landmask(lm_image, struct_elem; num_pixels_closing=num_pixels_closing)

        # test for percent difference in landmask images, ignore edges because we are not padding in Julia before applying strel_file
        @test (@test_approx_eq_sigma_eps masked_image[strel_h:end-strel_h, strel_w:end-strel_w] matlab_landmask[strel_h:end-50, strel_w:end-50] [0,0] 0.005) == nothing 
    end
end