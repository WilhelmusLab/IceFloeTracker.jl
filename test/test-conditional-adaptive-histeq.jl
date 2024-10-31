using IceFloeTracker:
    _get_false_color_cloudmasked,
    convert_to_255_matrix,
    adapthisteq,
    conditional_histeq,
    rgb2gray,
    to_uint8,
    imadjust

begin
    datadir = joinpath(@__DIR__, "test_inputs/")
    path_true_color_image = joinpath(
        datadir, "NE_Greenland_truecolor.2020162.aqua.250m.tiff"
    )
    path_false_color_image = joinpath(
        datadir, "NE_Greenland_reflectance.2020162.aqua.250m.tiff"
    )
    true_color_image = float64.(load(path_true_color_image))
    false_color_image = float64.(load(path_false_color_image))
    dilated_landmask = BitMatrix(load(joinpath(datadir, "matlab_landmask.png")))
end

function test_cloud_image_workflow()
    @testset "Prereq cloud image" begin
        false_color_cloudmasked = _get_false_color_cloudmasked(;
            false_color_image=false_color_image,
            prelim_threshold=110.0,
            band_7_threshold=200.0,
            band_2_threshold=190.0,
        )

        @test [sum(false_color_cloudmasked[i, :, :]) for i in 1:3] == [1_736_661_355, 5_997_708_807, 6_083_703_526]
    end
end

function test_adaphisteq()
    @testset "Adaptive histogram equalization" begin
        img = convert_to_255_matrix(testimage("cameraman"))
        img_eq = adapthisteq(img)
        @test sum(img_eq) == 32_387_397
    end
end

function test_conditional_adaptivehisteq()
    @testset "Conditional adaptivehisteq" begin
        clouds = _get_false_color_cloudmasked(;
            false_color_image=false_color_image,
            prelim_threshold=110.0,
            band_7_threshold=200.0,
            band_2_threshold=190.0,
        )

        clouds_red = clouds[1, :, :]
        clouds_red[dilated_landmask] .= 0

        @test sum(clouds_red) == 1_320_925_065

        # Using rblocks = 8, cblocks = 6
        true_color_eq = conditional_histeq(true_color_image, clouds_red, 8, 6)

        # This differs from MATLAB script due to disparity in the implementations
        # of the adaptive histogram equalization / diffusion functions
        # For the moment testing for regression
        @test sum(to_uint8(true_color_eq[:, :, 1])) == 6_372_159_606

        # Use custom tile size
        side_length = size(true_color_eq, 1) ÷ 8
        true_color_eq = conditional_histeq(true_color_image, clouds_red, side_length)
        @test sum(to_uint8(true_color_eq[:, :, 1])) == 6_328_796_398
    end
end

function test_rgb2gray()
    @testset "RGB to grayscale" begin
        g = rgb2gray(true_color_image)
        @test g[1000, 1000] == 92 && g[2000, 2000] == 206
    end
end

function test_imadjust()
    @testset "imadjust" begin
        Random.seed!(123)
        img = rand(0:255, 512, 512)
        @test sum(imadjust(img)) == 33_457_734
    end
end

@testset "Conditional adaptivehisteq" begin
    test_cloud_image_workflow()
    test_adaphisteq()
    test_conditional_adaptivehisteq()
    test_rgb2gray()
    test_imadjust()
end
