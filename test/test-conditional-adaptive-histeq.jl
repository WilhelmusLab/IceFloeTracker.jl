begin
    datadir = joinpath(@__DIR__, "test_inputs/")
    path_true_color_image = joinpath(datadir, "NE_Greenland_truecolor.2020162.aqua.250m.tiff")
    path_false_color_image = joinpath(datadir, "NE_Greenland_reflectance.2020162.aqua.250m.tiff")
    true_color_image = float64.(load(path_true_color_image))
    false_color_image = float64.(load(path_false_color_image))
    dilated_landmask = BitMatrix(load(joinpath(datadir, "matlab_landmask.png")))
end

function test_cloud_image_workflow()
    @testset "Prereq cloud image" begin
        redchannel_cahe = IceFloeTracker._get_red_channel_cloud_cae(
            false_color_image=false_color_image,
            landmask=dilated_landmask,
            prelim_threshold=110.0,
            band_7_threshold=200.0,
            band_2_threshold=190.0,
        )

        @test sum(redchannel_cahe) == 1_320_925_065
    end
end


function test_adaphisteq()
    @testset "Adaptive histogram equalization" begin
        img = IceFloeTracker.convert_to_255_matrix(testimage("cameraman"))
        img_eq = IceFloeTracker.adapthisteq(img)
        @test sum(img_eq) == 32_387_397
    end
end

function test_conditional_adaptivehisteq()
    @testset "Conditional adaptivehisteq" begin

        # Using rblocks = 8, cblocks = 6
        true_color_eq = IceFloeTracker.conditional_histeq(
            true_color_image,
            false_color_image,
            dilated_landmask,
            8,
            6)

        # This differs from MATLAB script due to disparity in the implementations
        # of the adaptive histogram equalization / diffusion functions
        # For the moment testing for regression
        @test sum(IceFloeTracker.to_uint8(true_color_eq[:, :, 1])) == 6_372_159_606


        # Use custom tile size
        side_length = size(true_color_eq, 1) ÷ 8
        true_color_eq = IceFloeTracker.conditional_histeq(
            true_color_image,
            false_color_image,
            dilated_landmask,
            side_length)
        @test sum(IceFloeTracker.to_uint8(true_color_eq[:, :, 1])) == 6_328_796_398
    end
end

@testset "Conditional adaptivehisteq" begin
    test_cloud_image_workflow()
    test_adaphisteq()
    test_conditional_adaptivehisteq()
end
